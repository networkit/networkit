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
    with open("api/DevGuide.rst", "w+") as fout:
        fout.write(".. _devGuide:\n\n" + devGuide)

def prepareBuild():
    insertNewsReferences()
    createDevGuide()

def cleanUp():
    removeNewsReferences()
