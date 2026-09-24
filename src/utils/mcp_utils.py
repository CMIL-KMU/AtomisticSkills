import sys
import os
import anyio


def setup_mcp_stdout():
    """
    Standardize MCP stdout/stderr redirection to prevent library noise
    from polluting the MCP JSON-RPC protocol.

    This function:
    1. Captures the original stdout (FD 1).
    2. Redirects system-level FD 1 to FD 2 (stderr) to catch C-level prints.
    3. Reassigns sys.stdout to sys.stderr at the Python level.
    4. Returns a handle to the original stdout for the MCP transport.
    """
    try:
        # Check if we are using the safe launcher which saved the pipe to a custom FD
        custom_fd_str = os.environ.get("MCP_STDOUT_FD")
        if custom_fd_str:
            try:
                mcp_stdout_fd = int(custom_fd_str)
                # Ensure this FD is open and writable
                # We do NOT need to dup2(2, 1) because the shell already did it!
                # But we might need to patch sys.stdout just in case Python reset it.
                sys.stdout = sys.stderr
                return os.fdopen(mcp_stdout_fd, "wb", buffering=0)
            except Exception as e:
                sys.stderr.write(
                    f"Warning: Failed to use MCP_STDOUT_FD={custom_fd_str}: {e}\n"
                )
                # Fallback to normal logic

        # 1. Save the REAL stdout (the one used for MCP communication)
        mcp_stdout_fd = os.dup(1)

        # 2. Redirect system-level FD 1 to stderr (FD 2)
        os.dup2(2, 1)

        # 3. Create a handle to the REAL pipe
        # We use a raw binary file object so we can re-wrap it correctly in the transport
        mcp_pipe_binary = os.fdopen(mcp_stdout_fd, "wb", buffering=0)

        # 4. Patch Python's sys.stdout to use stderr
        sys.stdout = sys.stderr

        return mcp_pipe_binary
    except Exception as e:
        sys.stderr.write(f"Warning: Failed to setup robust MCP stdout isolation: {e}\n")
        return None


def run_mcp_server(mcp, mcp_pipe_binary):
    """Hand the saved wire to SDK 2.2's public stdio runner.

    Imports run with both Python and native stdout redirected to stderr. Once
    startup is complete, the SDK claims fd 1 itself and isolates tool/child
    output for the lifetime of its transport. No private server API is used.
    """
    if mcp_pipe_binary is None:
        raise RuntimeError("MCP stdout isolation failed; refusing an unsafe transport")
    tool_lock = anyio.Lock()

    async def serialize_tools(context, call_next):
        # SDK 2 runs synchronous tools in threads. Model wrappers and the active
        # research directory are process-wide state: retain serial tool calls.
        if context.method == "tools/call":
            async with tool_lock:
                return await call_next(context)
        return await call_next(context)

    mcp.middleware.append(serialize_tools)
    sys.stdout.flush()
    previous_stdout = sys.stdout
    os.dup2(mcp_pipe_binary.fileno(), 1)
    sys.stdout = sys.__stdout__
    try:
        mcp.run(transport="stdio")
    finally:
        os.dup2(2, 1)
        sys.stdout = previous_stdout
        mcp_pipe_binary.close()
