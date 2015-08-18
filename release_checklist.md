# RELEASE CHECKLIST

## Development related

    1. Optional, if not already done during ongoing development:
        * Pull new features/code from forks into the Dev branch.
        * Merge branches into the Dev branch.
    2. Make sure, that the all the Unit tests (option '-t/--tests') run properly.
    3. If not done in the branches and forks: Document changes relevant to the user in the appropriate files (markdown, PDF, IPython notebook).
    4. Make sure, that the user guide Notebook runs properly.
    5. Update version number (setup.py, version.py).
    6. Merge Dev branch into default (release) branch.
    7. Set version tag in default branch.
    8. Merge default branch back into Dev branch.

## Website related

    * Write news item.
    * Update "getting started" page.
    * Update files (networkit.zip, PDF, user guide,...) and their links.
    * Update "documentation" page and C++ and Python documentations.

## Misc
    * Write a mail including release notes to the mailing list.
