"""Protect existing environments and retain failed installation evidence."""
from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

from tools import provision_environments as provision


def test_generated_mcp_config_uses_only_selected_environment(tmp_path):
    from configure_mcp import load_mcp_servers
    servers = load_mcp_servers(str(tmp_path))
    assert servers
    assert all(server["env"]["PYTHONNOUSERSITE"] == "1" for server in servers.values())
    assert all(server["command"].startswith(str(tmp_path)) for server in servers.values())


def test_existing_directory_is_never_modified(tmp_path, monkeypatch):
    item = provision.plan("base-agent")
    work = tmp_path / item["name"]
    work.mkdir()
    protected = work / "receipt.json"
    protected.write_text("previous evidence")
    monkeypatch.setattr(provision.subprocess, "run", lambda *a, **k: (_ for _ in ()).throw(AssertionError("Must not install")))
    assert provision.install(item, tmp_path, "conda")["status"] == "existing-prefix-or-receipt"
    assert protected.read_text() == "previous evidence"


def test_unavailable_local_source_blocks_before_conda(tmp_path, monkeypatch):
    item = deepcopy(provision.plan("base-agent"))
    item["missing_local_sources"] = ["unavailable provider checkout"]
    monkeypatch.setattr(provision.subprocess, "run", lambda *a, **k: (_ for _ in ()).throw(AssertionError("Must not install")))
    result = provision.install(item, tmp_path, "conda")
    assert result["status"] == "blocked-local-source"
    assert not Path(result["prefix"]).exists()


def test_failed_installer_retains_failure_and_source_identity(tmp_path, monkeypatch):
    item = provision.plan("base-agent")
    def fail(command, **kwargs):
        kwargs["stdout"].write("fixture dependency resolver failure\n")
        return SimpleNamespace(returncode=7)
    monkeypatch.setattr(provision.subprocess, "run", fail)
    result = provision.install(item, tmp_path, "conda")
    retained = json.loads((tmp_path/item["name"]/"receipt.json").read_text())
    assert retained == result and result["exit_code"] == 7
    assert result["status"] == "installation-failed"
    assert result["source_sha256"] == item["source_sha256"]
    assert "resolver failure" in (tmp_path/item["name"]/"install-0.log").read_text()
