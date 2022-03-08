import os
import pandas as pd
import tarfile
import gzip

def gunzip(source_filepath, dest_filepath, block_size=65536):
    with gzip.open(source_filepath, 'rb') as s_file, \
            open(dest_filepath, 'wb') as d_file:
        while True:
            block = s_file.read(block_size)
            if not block:
                break
            else:
                d_file.write(block)

test_path = "/data/Catherine/feh230"
extract_path = "/data/extract"
rewrite_path = "/data/rewrite"
for file in os.listdir(test_path):
    if file[-4:] == '.tar' and exists("{}/{}".format(extract_path,file[:-4]))==False and file[-9:-4] != 19925:
        my_tar = tarfile.open("{}/{}".format(test_path,file))
        my_tar.extractall("{}/{}".format(extract_path,file[:-4]))
        my_tar.close()
        os.mkdir("{}/{}".format(rewrite_path,file[:-4]))
        if len(os.listdir("{}/{}".format(extract_path,file[:-4]))) != 21 and len(os.listdir("{}/{}".format(extract_path,file[:-4]))) != 12:
            print(file)
        for gz in os.listdir("{}/{}".format(extract_path,file[:-4])):
            gunzip("{}/{}/{}".format(extract_path,file[:-4],gz),"{}/{}/{}".format(extract_path,file[:-4],gz[:-3]))
            os.remove("{}/{}/{}".format(extract_path,file[:-4],gz))
            f = open("{}/{}/{}".format(extract_path,file[:-4],gz[:-3]), "r")
            w = open("{}/{}/{}".format(rewrite_path,file[:-4],gz[:-3]), "w")
            lines_f = f.readlines()
            lines_w = lines_f
            if lines_w[2] == '\n':
                lines_w[3] = lines_f[3][:8] + '0' + lines_f[3][8:] 
                del lines_w[2]
            w.writelines(lines_w)
            f.close()
            w.close()
