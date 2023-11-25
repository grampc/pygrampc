import scikit_build_core.cmake as skbuild
from pathlib import Path

if __name__ == "__main__":
    # Create CMaker for the current directory with RELEASE build mode
    CMaker = skbuild.CMaker(skbuild.CMake.default_search(), Path.cwd(), Path.joinpath(Path.cwd(), "build"), "RELEASE")

    CMaker.configure()
    CMaker.build()
    print("Finished building")
