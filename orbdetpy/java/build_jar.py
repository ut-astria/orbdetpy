# build_jar.py - Rebuild Java sources and create JAR archive.
# Copyright (C) 2019 University of Texas
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from glob import iglob
from os import chdir, environ, getcwd, path, sep
from subprocess import run, PIPE, STDOUT

rootdir = path.dirname(path.dirname(path.abspath(__file__)))
srcdir = path.join(rootdir, "java")
libsdir = path.join(rootdir, "lib")
pushd = getcwd()
chdir(srcdir)

cpath = ""
csc = ":" if sep == "/" else ";"
for jar in iglob(path.join(libsdir, "**"), recursive = True):
    if (jar.endswith(".jar")):
        cpath += jar + csc
environ["CLASSPATH"] = cpath

for src in iglob("**", recursive = True):
    if (not src.endswith(".java")):
        continue

    print("Compiling " + src)
    res = run(["javac", "-Xlint:deprecation", "-Xlint:unchecked", src],
              stdout = PIPE, stderr = STDOUT, text = True)
    if (res.returncode != 0):
        print(res.stdout)
        chdir(pushd)
        exit(res.returncode)

cls = [f for f in iglob("**", recursive=True) if f.endswith(".class")]
jar = path.join(libsdir, "astria.jar")
print("Creating " + jar)
res = run(["jar","cf",jar,*cls], stdout=PIPE, stderr=STDOUT, text=True)
if (res.returncode != 0):
    print(res.stdout)
    chdir(pushd)
    exit(res.returncode)
chdir(pushd)
