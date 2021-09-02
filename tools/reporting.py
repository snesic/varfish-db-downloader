import time
import sys
import os
import csv
import hashlib

from awe import Page
from pytools.persistent_dict import PersistentDict

storage = PersistentDict("reportstorage")

now = time.time()
page = Page('VarFish DB Downloader')

FILES = []
ALL_DOWNLOAD_FILES_GRCH37 =[]
ALL_DOWNLOAD_FILES_GRCH38 = []
page.new("h1").new_text("Sanity Report")
print("Generating sanity reports start...", file=sys.stderr)

files = storage.fetch("files")
for file in files:
    if file.endswith('.md5'):
        continue
    print("path ..." + file, file=sys.stderr)

    if (file in FILES):
        continue
    FILES.append(file)
    print("file ..." + file, file=sys.stderr)

    md5File = file + '.md5'
    if not os.path.isfile(md5File):
        print("md5File file not found..." + md5File, file=sys.stderr)
        continue
    tsv_file = open(md5File)
    read_tsv = csv.reader(tsv_file, delimiter="\t")
    md5FromFile = ""
    for row in read_tsv:
        md5FromFile = row[0].split(" ")[0]
    tsv_file.close()

    with  open(file, "rb") as f:
        file_hash = hashlib.md5()
        while True:
            data = f.read(2 ** 20)
            if not data:
                break
            file_hash.update(data)
    digest = file_hash.hexdigest()
    equal = digest == md5FromFile

    if (file.startswith("GRCh37")):
        ALL_DOWNLOAD_FILES_GRCH37.append(['{}'.format(file), '{}'.format(digest),
                                          '{}'.format(md5FromFile), '{}'.format(equal)])
    else:
        ALL_DOWNLOAD_FILES_GRCH38.append(['{}'.format(file), '{}'.format(digest),
                                          '{}'.format(md5FromFile), '{}'.format(equal)])

grid = page.new_grid(columns=1, props={'gutter': 12})
grid.new_card('This is Sanity Checker description')
tabs = grid.new_tabs()
collapse = grid.new_collapse()
grid.new_divider()

tabs.new_tab('GrCh37').new_table(['File', 'Calculated Checksum', 'Checksum from MD5 file',
                                  'Result'], page_size=4).extend([
    ALL_DOWNLOAD_FILES_GRCH37[i]
    for i in range(len(ALL_DOWNLOAD_FILES_GRCH37))
])

tabs.new_tab('GrCh38').new_table(['File', 'Calculated Checksum', 'Checksum from MD5 file',
                                  'Result'], page_size=4).extend([
    ALL_DOWNLOAD_FILES_GRCH38[i]
    for i in range(len(ALL_DOWNLOAD_FILES_GRCH38))
])

page.start(True)

print("Generating sanity reports end...", file=sys.stderr)