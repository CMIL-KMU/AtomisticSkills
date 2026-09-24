"""Verify public SDK transport isolation and errors in fresh server processes."""
import asyncio
from pathlib import Path
import sys

import pytest


@pytest.mark.base
@pytest.mark.parametrize("mode", ["auto", "legacy"])
@pytest.mark.parametrize("launcher_fd", [False, True])
def test_native_and_python_output_never_pollute_mcp(tmp_path, mode, launcher_fd):
    from mcp import Client, StdioServerParameters

    root = Path(__file__).resolve().parents[2]
    server = tmp_path / "noisy_server.py"
    server.write_text(
        "import os, sys, time\n"
        f"sys.path.insert(0, {str(root)!r})\n"
        + ("os.environ['MCP_STDOUT_FD'] = str(os.dup(1))\nos.dup2(2, 1)\n" if launcher_fd else "")
        + "from src.utils.mcp_utils import setup_mcp_stdout, run_mcp_server\n"
        "pipe = setup_mcp_stdout()\n"
        "print('startup Python noise')\nos.write(1, b'startup native noise\\n')\n"
        "from mcp.server import MCPServer\n"
        "server = MCPServer('isolation-fixture')\n"
        "@server.tool()\ndef noisy(value: int) -> str:\n"
        "    print('tool Python noise')\n    os.write(1, b'tool native noise\\n')\n"
        "    if value < 0:\n        raise ValueError('fixture rejection')\n"
        "    return str(value)\n"
        "busy = False\n"
        "@server.tool()\ndef serial() -> str:\n"
        "    global busy\n"
        "    if busy:\n        raise RuntimeError('overlapping stateful calls')\n"
        "    busy = True\n    time.sleep(0.1)\n    busy = False\n    return 'done'\n"
        "run_mcp_server(server, pipe)\n"
    )

    async def exercise():
        params = StdioServerParameters(command=sys.executable, args=[str(server)])
        async with asyncio.timeout(20):
            async with Client(params, mode=mode) as client:
                listed = await client.list_tools()
                assert listed.tools[0].input_schema["properties"]["value"]["type"] == "integer"
                result = await client.call_tool("noisy", {"value": 7})
                assert not result.is_error and result.content[0].text == "7"
                error = await client.call_tool("noisy", {"value": -1})
                assert error.is_error
                assert error.model_dump(by_alias=True)["isError"] is True
                concurrent = await asyncio.gather(*(client.call_tool("serial", {}) for _ in range(3)))
                assert all(not r.is_error and r.content[0].text == "done" for r in concurrent)
    asyncio.run(exercise())
