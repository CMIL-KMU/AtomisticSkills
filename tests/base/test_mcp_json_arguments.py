"""Exercise the original server through MCP, including SDK JSON pre-parsing."""
import asyncio
from pathlib import Path
import sys

import pytest


@pytest.mark.base
@pytest.mark.parametrize("matrix", ['[2,1,1]', [2, 1, 1]])
def test_supercell_json_string_and_structured_matrix(tmp_path, matrix):
    from mcp import Client, StdioServerParameters
    from pymatgen.core import Lattice, Structure

    source = tmp_path / "input.cif"
    output = tmp_path / "supercell.cif"
    original = Structure(Lattice.cubic(4), ["Si"], [[0, 0, 0]])
    original.to(filename=str(source))
    root = Path(__file__).resolve().parents[2]

    async def exercise():
        server = StdioServerParameters(command=sys.executable,
            args=[str(root / "src/mcp_server/base_server.py")], cwd=str(tmp_path))
        async with asyncio.timeout(30):
            async with Client(server) as client:
                result = await client.call_tool("supercell_expansion", dict(
                        structure_path=str(source), scaling_matrix_json=matrix, save_to_file=str(output)))
                assert not result.is_error
                assert result.content[0].text.startswith("Successfully created supercell")
    asyncio.run(exercise())
    generated = Structure.from_file(output)
    assert len(generated) == 2
    assert generated.volume == pytest.approx(original.volume * 2)
