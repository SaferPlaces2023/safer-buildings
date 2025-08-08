# -------------------------------------------------------------------------------
# License:
# Copyright (c) 2012-2022 Luzzi Valerio
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#
#
# Name:        filesystem.py
# Purpose:
#
# Author:      Luzzi Valerio
#
# Created:     16/12/2019
# -------------------------------------------------------------------------------

import datetime
import hashlib
import os
import sys
import tempfile
from .module_log import Logger


_AUX_FILES = {
    'shp': ['shx', 'dbf', 'prj', 'cpg'],
}

_GPD_DRIVERS = lambda ext : {
    'geojson': 'GeoJSON',
    'json': 'GeoJSON',
    'shp': 'ESRI Shapefile',
    'gpkg': 'GPKG'
}.get(ext.lower())


def now():
    return datetime.datetime.now()


def total_seconds_from(t):
    return (datetime.datetime.now() - t).total_seconds()

def python_path():
    """
    python_path
    """
    pathname, _ = os.path.split(normpath(sys.executable))
    return pathname


def normpath(pathname):
    """
    normpath
    """
    if not pathname:
        return ""
    normp = os.path.normpath(pathname.replace("\\", "/")).replace("\\", "/")
    normp = f's3://{normp.lstrip("s3:").lstrip("/")}' if normp.startswith("s3:") else normp
    return normp


def juststem(pathname):
    """
    juststem
    """
    pathname = os.path.basename(pathname)
    root, _ = os.path.splitext(pathname)
    return root


def justpath(pathname, n=1):
    """
    justpath
    """
    for _ in range(n):
        pathname, _ = os.path.split(normpath(pathname))
    if pathname == "":
        return "."
    return normpath(pathname)


def justfname(pathname):
    """
    justfname - returns the basename
    """
    return normpath(os.path.basename(normpath(pathname)))


def justext(pathname):
    """
    justext
    """
    pathname = os.path.basename(normpath(pathname))
    _, ext = os.path.splitext(pathname)
    return ext.lstrip(".")


def forceext(pathname, newext):
    """
    forceext
    """
    root, _ = os.path.splitext(normpath(pathname))
    pathname = root + ("." + newext if len(newext.strip()) > 0 else "")
    return normpath(pathname)


def get_aux_files(pathname):
    """
    get_aux_files
    """
    pathname = normpath(pathname)
    ext = justext(pathname)
    aux_exts = _AUX_FILES.get(ext, [])
    return [forceext(pathname, aux_ext) for aux_ext in aux_exts]


def strtofile(text, filename, append=False):
    """
    strtofile
    """
    try:
        flag = "a" if append else "w"
        if isinstance(text, (str,)):
            text = text.encode("utf-8")
        if isinstance(text, (bytes,)):
            flag += 'b'
        mkdirs(justpath(filename))
        with open(filename, flag) as stream:
            if text:
                stream.write(text)
    except Exception as ex:
        Logger.error(ex)
        return ""
    return filename


def filetostr(filename):
    """
    filetostr
    """
    try:
        with open(filename, "r", encoding="utf-8") as stream:
            return stream.read()
    except:
        return None
    
    
def text_replace(filename, old, new):
    """
    text_replace
    """
    text = filetostr(filename)
    if text:
        text = text.replace(old, new)
        strtofile(text, filename)
        return True
    return False


def isfile(pathname):
    """
    isfile
    """
    return pathname and isinstance(pathname, str) and os.path.isfile(pathname)


def israster(pathname):
    """
    israster
    """
    return isfile(pathname) and pathname.lower().endswith(".tif")


def isshape(pathname):
    """
    isshape
    """
    return isfile(pathname) and pathname.lower().endswith(".shp")


def iszip(pathname):
    """
    iszip
    """
    if isfile(pathname):
        with open(pathname, 'rb') as file:
            # Read the first 4 bytes
            header = file.read(4)
            # Check if the file header matches the ZIP file signature
            return  header == b'PK\x03\x04'
    return False
    
def remove(pathname):
    """
    remove
    """
    try:
        if isshape(pathname):
            for ext in ["shp", "shx", "dbf", "prj", "cpg", "mta"]:
                if isfile(forceext(pathname, ext)):
                    os.unlink(forceext(pathname, ext))
        elif israster(pathname):
            for ext in ["tif", "tfw", "aux.xml"]:
                if isfile(forceext(pathname, ext)):
                    os.unlink(forceext(pathname, ext))
        elif isfile(pathname):
            os.unlink(pathname)
        
        return True
    except Exception as ex:
        Logger.error(ex)
        return False
    

def mkdirs(pathname):
    """
    mkdirs - create a folder
    """
    try:
        if os.path.isfile(pathname):
            pathname = justpath(pathname)
        os.makedirs(pathname)
    except:
        pass
    return os.path.isdir(pathname)


def tempdir(suffix=""):
    """
    tempdir
    """
    res = f"{tempfile.gettempdir()}/{suffix}"
    os.makedirs(res, exist_ok=True)
    return res


def tempfilename(prefix="", suffix=""):
    """
    return a temporary filename
    """
    return normpath(tempfile.gettempdir() + "/" + datetime.datetime.strftime(now(), f"{prefix}%Y%m%d%H%M%S%f{suffix}"))


def listify(text, sep=",", trim=False):
    """
    listify -  make a list from string
    """
    if text is None:
        return []
    elif isinstance(text, str):
        arr = text.split(sep)
        if trim:
            arr = [item.strip() for item in arr]
        return arr
    elif isinstance(text, (tuple, list)):
        return text
    return [text]


def md5sum(filename):
    """
    md5sum - returns themd5 of the file
    """
    if os.path.isfile(filename):
        f = open(filename, mode='rb')
        d = hashlib.md5()
        while True:
            buf = f.read(4096)
            if not buf:
                break
            d.update(buf)
        f.close()
        res = d.hexdigest()
        return res
    else:
        return ""


def md5text(text):
    """
    md5text - Returns the md5 of the text
    """
    if text!=None:
        hash = hashlib.md5()
        if isinstance(text, (bytes, bytearray)):
            hash.update(text)
        else:
            hash.update(text.encode("utf-8"))
        return hash.hexdigest()
    return None

