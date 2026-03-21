import json
import shutil

from pathlib import Path


def main():
    names = [
        "BCC-Fe-Fe",
        "BCC-Fe-Mn",
        "BCC-Fe-Ni",
        "BCC-Mn-Mn",
        "BCC-Mn-Ni",
        "BCC-Ni-Ni",
        "FCC-Fe-Fe",
        "FCC-Fe-Mn",
        "FCC-Fe-Ni",
        "FCC-Mn-Fe",
        "FCC-Mn-Mn",
        "FCC-Mn-Ni",
        "FCC-Ni-Fe",
        "FCC-Ni-Mn",
        "FCC-Ni-Ni",
        "HCP-Fe-Fe",
        "HCP-Fe-Mn",
        "HCP-Fe-Ni",
        "HCP-Mn-Fe",
        "HCP-Mn-Mn",
        "HCP-Mn-Ni",
        "HCP-Ni-Fe",
        "HCP-Ni-Mn",
        "HCP-Ni-Ni",
        "SER-Fe",
        "SER-Mn",
        "SER-Ni",
    ]

    for name in names:
        workdir = Path("data/EndMembers").joinpath(name)
        json_path = workdir.joinpath(f"{name}-qha.json")

        if not workdir.exists():
            continue

        # Skip if already done
        if json_path.exists():
            with open(json_path, "r", encoding="utf-8") as jf:
                result = json.load(jf)
            if result.get("state") == "successful":
                # print(f"System {name} already done")
                continue

        print(f"System {name} failed")
        shutil.rmtree(workdir)


if __name__ == "__main__":
    main()
