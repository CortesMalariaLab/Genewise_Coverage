import os

def list_dir(in_dir):
    try:
        list_files = os.listdir(in_dir)
        return list_files
    except OSError as e:
        return []

def create_dir(out_dir):
    try:
        os.makedirs(out_dir)
    except OSError as e:
        return False

def iter_files(lf, in_dir, out_dir, thr):
    create_dir(out_dir)
    for file in lf:
        old = in_dir +'/'+file
        new = out_dir+'/'+file

        new_handle = open(new, 'w')
        with open(old, 'r') as old_handle:
            for line in old_handle:
                nl = new_line(line, thr)
                new_handle.write(nl)
            new_handle.close()

def new_line(line, thr):
    lls = line.rstrip().split('\t')
    if float(lls[3]) < thr:
        lls[3] = thr
    return f'{lls[0]}\t{lls[1]}\t{lls[2]}\t{lls[3]}\n'

def main_function(in_dir, out_dir, thr):
    ld = list_dir(in_dir)
    iter_files(ld, in_dir, out_dir, thr)

if __name__ == '__main__':
    ## Input Dir (Tracks Normalitzats per input)
    in_dir = '/PROJECTES/MALARIA_EPIGENETICS/Projects/ReCHIP_GDV1_Project_ene25/ALL_Tracks/NormInput/'

    ##Output Dir
    out_dir = '/PROJECTES/MALARIA_EPIGENETICS/Projects/ReCHIP_GDV1_Project_ene25/NEWGenewise_Coverage/ReNormTracks/'

    ##Threshold
    thr = -1.0

    main_function(in_dir, out_dir, thr)

