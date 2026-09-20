#!/usr/bin/env python3
"""Serial, package-only orchestration; no mathematical predicates live here.

Python 3.8+ standard library. See README.md for the root manifest contract.
"""

import argparse
import datetime
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import signal
import stat
import string
import subprocess
import sys
import tempfile
import time
import uuid


SCHEMA = "fujita-supplement-root-v1"
REPORT_SCHEMA = "fujita-supplement-run-v1"
MAX_JSON_BYTES = 64 * 1024 * 1024
BOUNDARY = (
    "Documentary file/source/version checks and execution of declared section "
    "validators only. No formal proof, geometric applicability verdict, native "
    "paper-verifier acceptance, or certification of the historical fact graph."
)


class Invalid(Exception):
    pass


def require(condition, message):
    if not condition:
        raise Invalid(message)


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def nonempty_string(value):
    return isinstance(value, str) and bool(value.strip())


def strict_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, "duplicate JSON key: " + key)
        result[key] = value
    return result


def reject_constant(value):
    raise Invalid("non-finite JSON constant: " + value)


def finite_float(value):
    result = float(value)
    require(math.isfinite(result), "non-finite JSON number: " + value)
    return result


def read_json(path):
    require(path.is_file(), "missing JSON file: " + str(path))
    require(path.stat().st_size <= MAX_JSON_BYTES, "JSON file exceeds 64 MiB: " + str(path))
    with path.open(encoding="utf-8") as stream:
        return json.load(stream, object_pairs_hook=strict_object, parse_constant=reject_constant,
                         parse_float=finite_float)


def identity(path):
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
            size += len(block)
    return {"sha256": digest.hexdigest(), "bytes": size}


def relative(value):
    require(nonempty_string(value) and "\\" not in value, "expected a portable relative path")
    path = PurePosixPath(value)
    require(not path.is_absolute() and all(part not in ("", ".", "..") for part in value.split("/")),
            "noncanonical or escaping relative path: " + value)
    require(not re.match(r"^[A-Za-z]:", value), "drive-qualified path: " + value)
    return path


def local_path(base, value):
    rel = relative(value)
    path = base
    for part in rel.parts:
        path = path / part
        require(not path.is_symlink(), "symlink is not a portable package member: " + str(path))
    return path


def below(path, directory):
    return path == directory or directory in path.parents


def inventory(directory):
    require(directory.is_dir(), "missing section directory: " + str(directory))
    result = {}
    for current, directories, files in os.walk(directory, followlinks=False):
        for name in directories:
            require(not (Path(current) / name).is_symlink(), "symlink directory in section: " + name)
        for name in files:
            path = Path(current) / name
            require(stat.S_ISREG(path.lstat().st_mode), "nonregular section member: " + str(path))
            result[path.relative_to(directory).as_posix()] = identity(path)
    return dict(sorted(result.items()))


def json_pointer(document, pointer):
    require(isinstance(pointer, str) and (pointer == "" or pointer.startswith("/")),
            "invalid JSON pointer")
    node = document
    for token in pointer.split("/")[1:] if pointer else []:
        require(not re.search(r"~(?![01])", token), "invalid JSON pointer escape")
        token = token.replace("~1", "/").replace("~0", "~")
        if isinstance(node, list):
            require(re.fullmatch(r"0|[1-9][0-9]*", token) is not None, "invalid array pointer: " + pointer)
            index = int(token)
            require(index < len(node), "missing JSON pointer: " + pointer)
            node = node[index]
        else:
            require(isinstance(node, dict) and token in node, "missing JSON pointer: " + pointer)
            node = node[token]
    return node


def assertions(document, checks, mandatory=False):
    require(isinstance(checks, list) and (not mandatory or checks), "nonempty report_checks required")
    for check in checks:
        require(isinstance(check, dict) and "pointer" in check and "equals" in check,
                "an assertion needs pointer and equals")
        observed = json_pointer(document, check["pointer"])
        # Structural assertions let a report retain variable timing metadata
        # while still binding its exact family count and complete field set.
        if "measure" in check:
            measure = check["measure"]
            require(measure in ("array_length", "object_keys"), "unsupported assertion measure")
            if measure == "array_length":
                require(isinstance(observed, list), "array_length needs an array")
                observed = len(observed)
            else:
                require(isinstance(observed, dict), "object_keys needs an object")
                observed = sorted(observed)
        # JSON types matter: True is not the integer 1.
        require(json.dumps(observed, sort_keys=True) == json.dumps(check["equals"], sort_keys=True),
                "JSON assertion failed at " + check["pointer"])


def explicit_failures(document):
    """Catch contradictory conventional failure signals as well as pinned checks."""
    failures = []
    stack = [("", document)]
    bad_status = {"fail", "failed", "failure", "error", "invalid", "rejected", "aborted",
                  "timeout", "incomplete", "skipped"}
    while stack:
        pointer, node = stack.pop()
        if isinstance(node, dict):
            for key, value in node.items():
                location = pointer + "/" + key.replace("~", "~0").replace("/", "~1")
                folded = key.lower()
                bad = folded in {"ok", "success", "passed"} and value is False
                bad |= folded in {"status", "verdict"} and isinstance(value, str) and value.lower() in bad_status
                bad |= folded in {"error", "errors", "failure", "failures", "failed"} and bool(value)
                if bad:
                    failures.append(location)
                stack.append((location, value))
        elif isinstance(node, list):
            stack.extend((pointer + "/" + str(index), value) for index, value in enumerate(node))
    return sorted(failures)


def report_failures(document, entry):
    """Permit explicitly pinned mathematical error allowances, never diagnostics.

    Some native reports call their rational error budget `error`. A root
    binding may classify only exact nonnegative rational values at exact
    pointers, with the same assertion also required in report_checks.
    Every other conventional failure remains a failure.
    """
    allowances = entry.get("report_error_allowances", [])
    require(isinstance(allowances, list), "report_error_allowances must be a list")
    ignored = set()
    for allowance in allowances:
        require(isinstance(allowance, dict) and set(allowance) == {"pointer", "equals"},
                "error allowance needs exactly pointer and equals")
        pointer, value = allowance["pointer"], allowance["equals"]
        require(isinstance(pointer, str) and pointer.endswith("/error") and pointer not in ignored,
                "error allowance requires a unique exact /error pointer")
        require(isinstance(value, str) and re.fullmatch(r"(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?", value),
                "error allowance must be an explicit nonnegative rational string")
        require(allowance in entry.get("report_checks", []),
                "error allowance must also have an exact report assertion")
        assertions(document, [allowance])
        ignored.add(pointer)
    return [pointer for pointer in explicit_failures(document) if pointer not in ignored]


def declaration(manifest, requested):
    require(isinstance(manifest, dict), "root manifest must be a JSON object")
    require(manifest.get("schema") == SCHEMA, "unsupported root manifest schema")
    require(manifest.get("ready") is True, "root manifest is a draft: integration is not complete")
    require(nonempty_string(manifest.get("package_version")), "missing package_version")
    expected = manifest.get("expected_sections")
    require(isinstance(expected, list) and expected, "expected_sections must not be empty")
    require(all(isinstance(s, str) and re.fullmatch(r"[a-zA-Z0-9][a-zA-Z0-9_-]*", s) for s in expected),
            "invalid section identifier")
    require(len(expected) == len(set(expected)), "duplicate expected section")
    entries = manifest.get("sections")
    require(isinstance(entries, list) and entries, "missing section declarations")
    require(all(isinstance(entry, dict) and entry.get("id") in expected for entry in entries),
            "undeclared or malformed section entry")
    indexed = {entry["id"]: entry for entry in entries}
    require(len(indexed) == len(entries), "duplicate section declaration")
    require(set(indexed) == set(expected), "expected sections and declarations do not match")
    require(not requested or len(requested) == len(set(requested)), "duplicate --section selection")
    require(not requested or all(s in indexed for s in requested), "unknown or empty --section selection")
    selected = [s for s in expected if not requested or s in requested]
    require(selected, "section selection must not be empty")
    paths = [relative(indexed[s].get("path")) for s in expected]
    for index, path in enumerate(paths):
        require(not any(below(path, other) or below(other, path) for other in paths[:index]),
                "overlapping section directories")
    return expected, indexed, selected


def preflight(base, entry):
    directory = local_path(base, entry["path"])
    require(nonempty_string(entry.get("version")), "missing section version")
    manifest_path = local_path(directory, entry.get("manifest"))
    checker = local_path(directory, entry.get("checker"))
    data = local_path(directory, entry.get("data"))
    require(checker.is_file(), "missing checker")
    require(local_path(directory, entry.get("readme", "README.md")).is_file(), "missing section README")
    require(data.is_dir(), "missing data directory")
    section_manifest = read_json(manifest_path)
    require(isinstance(section_manifest, dict) and section_manifest, "empty section manifest")
    assertions(section_manifest, entry.get("manifest_checks", []))
    declared = entry.get("files")
    require(isinstance(declared, list) and declared, "missing file identities")
    pinned = {}
    for member in declared:
        require(isinstance(member, dict), "malformed file identity")
        name = relative(member.get("path")).as_posix()
        require(name not in pinned, "duplicate file identity: " + name)
        require(isinstance(member.get("sha256"), str) and re.fullmatch(r"[0-9a-f]{64}", member["sha256"]),
                "missing/invalid SHA256: " + name)
        require(type(member.get("bytes")) is int and member["bytes"] >= 0, "invalid byte count: " + name)
        pinned[name] = {"sha256": member["sha256"], "bytes": member["bytes"]}
    observed = inventory(directory)
    require(set(observed) == set(pinned), "file inventory mismatch; missing=" + str(sorted(set(pinned) - set(observed)))
            + "; undeclared=" + str(sorted(set(observed) - set(pinned))))
    require(observed == pinned, "declared file identity mismatch")
    data_members = {name for name in observed if below(PurePosixPath(name), PurePosixPath(entry["data"]))}
    require(data_members, "empty data directory")
    mappings = entry.get("source_mappings")
    require(isinstance(mappings, list) and mappings, "missing source mappings")
    families, mapped = set(), set()
    for mapping in mappings:
        require(isinstance(mapping, dict), "malformed source mapping")
        family = mapping.get("family")
        require(nonempty_string(family) and family.startswith(entry["id"] + "/")
                and len(family) > len(entry["id"]) + 1, "family must be section-id/stable-key")
        require(family not in families, "duplicate source-mapping family")
        families.add(family)
        require(nonempty_string(mapping.get("source_version")), "missing mapped source version")
        references = mapping.get("source_refs")
        require(isinstance(references, list) and references and all(nonempty_string(ref) for ref in references),
                "missing documentary source references")
        payloads = mapping.get("data_files")
        require(isinstance(payloads, list) and payloads and all(isinstance(name, str) for name in payloads),
                "missing mapped data files")
        require(len(payloads) == len(set(payloads)) and set(payloads) <= data_members,
                "duplicate or nonexistent mapped data file")
        mapped.update(payloads)
    require(mapped == data_members, "source mappings do not cover all delivered data files")
    checks = entry.get("report_checks")
    require(isinstance(checks, list) and checks, "explicit report success assertions are required")
    for check in checks:
        require(isinstance(check, dict) and isinstance(check.get("pointer"), str) and "equals" in check,
                "malformed report success assertion")
    command = entry.get("command")
    require(isinstance(command, list) and len(command) >= 4 and all(nonempty_string(arg) for arg in command),
            "explicit command argument array required")
    require(command[:2] == ["{python}", "{checker}"], "command must start with {python}, {checker}")
    fields = []
    for argument in command:
        for _, field, spec, conversion in string.Formatter().parse(argument):
            if field is not None:
                require(field in {"python", "checker", "data", "report", "manifest"}
                        and not spec and conversion is None, "unsupported command placeholder")
                fields.append(field)
    require(fields.count("data") == 1 and fields.count("report") == 1,
            "command needs exactly one data and one fresh report argument")
    timeout = entry.get("timeout_seconds", 3600)
    require(type(timeout) in (float, int) and math.isfinite(timeout) and timeout > 0,
            "timeout_seconds must be finite and positive")
    return directory, observed, timeout


def stop_process(process):
    if os.name == "posix":
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    else:
        process.kill()
    process.wait()


def run_section(base, entry, directory, pinned, timeout, run_directory):
    output = run_directory / entry["id"]
    output.mkdir()
    report_path = output / "report.json"
    replacements = {"python": sys.executable, "checker": str(local_path(directory, entry["checker"])),
                    "data": str(local_path(directory, entry["data"])), "report": str(report_path),
                    "manifest": str(local_path(directory, entry["manifest"]))}
    command = [arg.format(**replacements) for arg in entry["command"]]
    # Disable bytecode output and inherited Python options (including -O).
    command[1:1] = ["-B", "-E"]
    result = {"id": entry["id"], "status": "fail", "ran": False, "version": entry["version"],
              "source_mappings": entry["source_mappings"], "documentary_prerequisites": "pass",
              "input_files": pinned, "command": command, "cwd": str(directory),
              "timeout_seconds": timeout, "report_path": str(report_path), "errors": []}
    errors = result["errors"]
    started_ns = time.time_ns()
    result["started_at"] = utc()
    start = time.monotonic()
    try:
        require(inventory(directory) == pinned, "section inputs changed after preflight")
        require(not report_path.exists(), "report path must be absent before launch")
        with (output / "stdout.log").open("wb") as stdout, (output / "stderr.log").open("wb") as stderr:
            process = subprocess.Popen(command, cwd=str(directory), stdin=subprocess.DEVNULL,
                                       stdout=stdout, stderr=stderr, start_new_session=(os.name == "posix"))
            result["ran"] = True
            try:
                result["exit_code"] = process.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                stop_process(process)
                result["exit_code"] = process.returncode
                errors.append("section timeout")
            except BaseException:
                stop_process(process)
                raise
        if result["exit_code"] != 0:
            errors.append("checker exited nonzero: " + str(result["exit_code"]))
        require(report_path.exists(), "checker did not create its fresh report")
        info = report_path.lstat()
        require(stat.S_ISREG(info.st_mode) and info.st_nlink == 1, "report must be a regular, unlinked file")
        # Two seconds allows coarse filesystem timestamp resolution. The new,
        # unpredictable directory is the primary freshness guarantee.
        require(info.st_mtime_ns >= started_ns - 2_000_000_000, "stale report modification time")
        result["report_identity"] = identity(report_path)
        payload = read_json(report_path)
        require(isinstance(payload, dict) and payload, "report must be a nonempty JSON object")
        result["section_report"] = payload
        failures = report_failures(payload, entry)
        require(not failures, "report contains failure signals: " + str(failures))
        assertions(payload, entry["report_checks"], mandatory=True)
    except KeyboardInterrupt:
        errors.append("interrupted during section execution")
        result["interrupted"] = True
        if result["ran"]:
            result["exit_code"] = process.returncode
    except (Invalid, OSError, ValueError, TypeError, RecursionError, OverflowError) as error:
        errors.append(str(error))
    finally:
        result["finished_at"] = utc()
        result["elapsed_seconds"] = round(time.monotonic() - start, 6)
        try:
            after = inventory(directory)
            result["inputs_unchanged"] = after == pinned
            if after != pinned:
                result["input_files_after"] = after
                errors.append("checker inputs changed during execution")
        except (Invalid, OSError) as error:
            errors.append("post-run identity check: " + str(error))
        result["logs"] = {name: {"path": str(output / name), **identity(output / name)}
                          for name in ("stdout.log", "stderr.log") if (output / name).is_file()}
    result["status"] = "pass" if not errors else "fail"
    return result


def write_report(path, document):
    path.parent.mkdir(parents=True, exist_ok=True)
    # Refuse to overwrite old aggregates as well as input files. A caller must
    # choose a new destination; previous evidence is never silently replaced.
    with path.open("x", encoding="utf-8") as stream:
        json.dump(document, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default="manifest.json", help="root JSON manifest; paths resolve beside it")
    parser.add_argument("--section", action="append", help="exact section ID; repeat to select several")
    parser.add_argument("--report", "--output", dest="report", default="root-report.json",
                        help="new aggregate JSON file, relative to the current directory")
    args = parser.parse_args(argv)
    aggregate = {"schema": REPORT_SCHEMA, "run_id": str(uuid.uuid4()), "status": "fail",
                 "started_at": utc(), "claim_boundary": BOUNDARY, "sections": [], "errors": [],
                 "selected_sections": [], "expected_sections": [], "ran_sections": [],
                 "full_package_run": False, "python": {"executable": sys.executable, "version": sys.version}}
    manifest_path = Path(args.manifest).resolve()
    report_path = Path(args.report).absolute()
    can_write = False
    exit_code = 2
    try:
        require(not report_path.exists() and not report_path.is_symlink(), "aggregate report already exists")
        require(report_path.resolve() != manifest_path, "aggregate report collides with root manifest")
        can_write = True
        manifest = read_json(manifest_path)
        aggregate["root_manifest"] = {"path": str(manifest_path), **identity(manifest_path)}
        aggregate["runner"] = identity(Path(__file__))
        # Record expected IDs even for a deliberate draft rejection.
        if isinstance(manifest, dict):
            aggregate["expected_sections"] = manifest.get("expected_sections", [])
            aggregate["package_version"] = manifest.get("package_version")
            for entry in manifest.get("sections", []) if isinstance(manifest.get("sections"), list) else []:
                if isinstance(entry, dict) and nonempty_string(entry.get("path")):
                    section_directory = local_path(manifest_path.parent, entry["path"]).resolve()
                    if below(report_path.resolve(), section_directory):
                        can_write = False
                        raise Invalid("aggregate/output directory must be outside section inputs")
        expected, entries, selected = declaration(manifest, args.section)
        aggregate["selected_sections"] = selected
        aggregate["full_package_run"] = selected == expected
        report_path.parent.mkdir(parents=True, exist_ok=True)
        run_directory = Path(tempfile.mkdtemp(prefix="section-runs-", dir=str(report_path.parent))).resolve()
        aggregate["run_directory"] = str(run_directory)
        exit_code = 0
        for section_id in expected:
            if section_id not in selected:
                aggregate["sections"].append({"id": section_id, "status": "not_selected", "ran": False})
                continue
            entry = entries[section_id]
            try:
                directory, pinned, timeout = preflight(manifest_path.parent, entry)
                result = run_section(manifest_path.parent, entry, directory, pinned, timeout, run_directory)
            except (Invalid, OSError, ValueError, TypeError, RecursionError, OverflowError) as error:
                result = {"id": section_id, "status": "fail", "ran": False,
                          "documentary_prerequisites": "fail", "errors": [str(error)]}
            aggregate["sections"].append(result)
            if result["ran"]:
                aggregate["ran_sections"].append(section_id)
            if result["status"] != "pass":
                exit_code = 1
            print(section_id + ": " + result["status"], flush=True)
            if result.get("interrupted"):
                aggregate["errors"].append("interrupted; remaining sections were not run")
                exit_code = 130
                processed = {row["id"] for row in aggregate["sections"]}
                aggregate["sections"].extend({"id": sid, "status": "not_run" if sid in selected else "not_selected",
                                              "ran": False} for sid in expected if sid not in processed)
                break
        require(identity(manifest_path) == {key: aggregate["root_manifest"][key] for key in ("sha256", "bytes")},
                "root manifest changed during execution")
        aggregate["status"] = "pass" if exit_code == 0 else "fail"
    except KeyboardInterrupt:
        aggregate["errors"].append("interrupted; aggregate is incomplete")
        exit_code = 130
    except (Invalid, OSError, ValueError, TypeError, RecursionError, OverflowError) as error:
        aggregate["errors"].append(str(error))
        exit_code = 2
    aggregate["exit_code"] = exit_code
    aggregate["finished_at"] = utc()
    if can_write:
        try:
            write_report(report_path, aggregate)
        except (OSError, ValueError) as error:
            print("cannot write aggregate: " + str(error), file=sys.stderr)
            return 2
    if aggregate["errors"]:
        print("; ".join(aggregate["errors"]), file=sys.stderr)
    print("aggregate: " + aggregate["status"] + "; report=" + str(report_path), flush=True)
    return exit_code


def terminate_as_interrupt(_signum, _frame):
    # Resource wrappers terminate the root with SIGTERM. Route that through
    # the same child-group cleanup and incomplete-report path as Ctrl-C.
    raise KeyboardInterrupt


if __name__ == "__main__":
    signal.signal(signal.SIGTERM, terminate_as_interrupt)
    sys.exit(main())
