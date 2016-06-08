# script that inserts news references into news.rst for sphinx
import re
import os
from subprocess import call

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
    insertPublicationsTags()

def cleanUp():
    removeNewsReferences()
    removePublicationsTags()
