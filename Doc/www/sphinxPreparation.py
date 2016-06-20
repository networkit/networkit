# script that prepares some files for parsing with sphinx
import re
import os
from subprocess import call
from shutil import copyfile

def insertNewsReferences():
    newsIdx = 1
    with open("news.rst") as fin, open("temp_news.rst", "w") as fout:
        currentLine = next(fin)
        for line in fin:
            if re.match('[~]+$', line):
                fout.write(".. _news-" + str(newsIdx) + ":\n\n")
                newsIdx = newsIdx+1
            fout.write(currentLine)
            currentLine = line
        fout.write(currentLine)
    os.rename("temp_news.rst", "news.rst")

def removeNewsReferences():
    with open("news.rst") as fin, open("temp_news.rst", "w") as fout:
        for line in fin:
            if not line.startswith(".. _news"):
                fout.write(line)
            else:
                next(fin)
    os.rename("temp_news.rst", "news.rst")

def createDevGuide():
    call(["pandoc", "--from=markdown",  "--to=rst", "--output=api/DevGuide.rst", "../DevGuide.mdown"])
    with open("api/DevGuide.rst", "r") as fin:
        devGuide = fin.read()

    devGuide = devGuide.splitlines()
    with open("api/DevGuide.rst", "w") as fout:
        fout.write(".. _devGuide:\n\n")
        for line in devGuide:
            if line == "Unit Tests and Testing":
                fout.write(".. _devGuide-unitTests:\n\n")
            fout.write(line + "\n")

def createReadme():
    # create the README.rst from the get_started.rst
    copyfile("get_started.rst", "../../README.rst")

    with open("../../README.rst", "r") as fin:
        readme = fin.read()

    readme = readme.splitlines()
    with open("../../README.rst", "w") as fout:
        regex = re.compile("<div")
        for line in readme:
            if not "|separator|" in line and regex.search(line) == None:
                fout.write(line + "\n")

    # Create the ../Readme.pdf
    #
    # pandoc can't parse internal references from .rst files so we convert
    # the .rst file to a .mdown file for creating the pdf file.
    with open("get_started.rst", "r") as fin:
        readme = fin.read()

    readme = readme.splitlines()
    with open("Readme.mdown", "w") as fout:
        regex = re.compile("`[^<>`]+`_", re.IGNORECASE)
        for line in readme:
            match = regex.search(line, re.IGNORECASE)
            while match:
                link = match.group(0)[1:-2]
                for i in range(0, len(readme)):
                    if readme[i].startswith(".. _" + link + ":"):
                        for j in range(i+1, len(readme)):
                            if readme[j]:
                                replStr = readme[j].lower()
                                replStr = replStr.replace(" ", "-")
                                line = regex.sub("`" + link + " <#" + replStr + ">`_", line, 1)
                                break
                        break
                    elif readme[i].startswith(link):
                        replStr = readme[i].lower()
                        replStr = replStr.replace(" ", "-")
                        line = regex.sub("`" + link + " <#" + replStr + ">`_", line, 1)
                        break
                match = regex.search(line, re.IGNORECASE)
            fout.write(line + "\n")

    call(["pandoc", "--from=rst",  "--to=markdown", "--output=Readme.mdown", "Readme.mdown"])
    with open("Readme.mdown", "r") as fin:
        readme = fin.read()

    readme = readme.splitlines()
    with open("Readme.mdown", "w") as fout:
        regex = re.compile("{.sourceCode}", re.IGNORECASE)
        for line in readme:
            if not "|separator|" in line:
                line = regex.sub("", line)
                fout.write(line + "\n")

    call(["pandoc", "-Vgeometry:margin=1.75cm", "--variable", "urlcolor=cyan", "Readme.mdown", "-o../Readme.pdf"])
    os.remove("Readme.mdown")


def insertPublicationsTags():
    tagIdx = 1
    increment = False
    with open("publications.rst") as fin:
        publications = fin.read()

    publications = publications.splitlines()
    with open("publications.rst", "w") as fout:
        for line in publications:
            while "\"collapseDiv\"" in line or "\"#collapseDiv\"" in line:
                line = line.replace("collapseDiv", "collapseDiv" + str(tagIdx), 1)
                if increment:
                    tagIdx += 1
                    increment = False
                else:
                    increment = True
            fout.write(line + "\n")

def removePublicationsTags():
    with open("publications.rst") as fin:
        publications = fin.read()

    publications = publications.splitlines()
    with open("publications.rst", "w") as fout:
        regex = re.compile("collapseDiv[\d]+", re.IGNORECASE)
        for line in publications:
            line = regex.sub("collapseDiv", line)
            fout.write(line + "\n")





def prepareBuild():
    insertNewsReferences()
    createDevGuide()
    createReadme()
    insertPublicationsTags()

def cleanUp():
    removeNewsReferences()
    removePublicationsTags()
