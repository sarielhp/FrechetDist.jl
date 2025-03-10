#!/usr/bin/python3

import csv
import os

deleted_trajectories = 0

def dropHeader(reader):
    for i in range(0, 6):
        next(reader)

def nearBeijing(lat, lon):
    return lat >= 39. and lat <= 41. and lon >= 115. and lon <= 117.

def parseAndWriteContent(reader, new_filename, filelog, filename):
    # create file content
    rows = []
    f_near = False
    for row in reader:
        rows.append([row[0], row[1]])
        # only use if near Beijing
        if nearBeijing(float(row[0]), float(row[1])):
            f_near = True;
        else:
            fnear = False

    global deleted_trajectories
    if  not f_near:
        new_filename = os.path.splitext( new_filename )[0]+ "_" + str( deleted_trajectories ) + ".txt"
        #+'.jpg'  # /home/user/somefile.jpg
        #new_filename = new_filename + "_" + str( deleted_trajectories )
        deleted_trajectories += 1

    print( filename, " ===> ", new_filename, "\n", file=filelog );
    # write file content
    with open(new_filename, 'w') as f:
        writer = csv.writer(f, delimiter=' ', quotechar='|')
        writer.writerows(rows)

        f.close()

    return f_near

def convert(filename, new_filename, filelog):
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter=',', quotechar='|')
        dropHeader(reader)
        success = parseAndWriteContent(reader, new_filename, filelog, filename)

        f.close()
        return success

def main():
    print("Start converting the trajectories. This takes some time...")

    filelog = open( 'conversion.log', 'w')
    all_files = []
    for subdir, dirs, files in os.walk('Data/'):
        for file in files:
            filename = subdir + os.sep + file
            if filename.endswith(".plt"):
                all_files.append(filename)

    all_files.sort()

    file_number = 0
    for filename in all_files:
        new_filename = 'data/' + str(file_number) + ".txt"
        success = convert(filename, new_filename, filelog)
        if success:
            file_number += 1

    filelog.close();
    print("Finished conversion of data.")

    global deleted_trajectories
    print("Deleted " + str(deleted_trajectories) + " trajectories.")

if __name__ == "__main__":
    main()
