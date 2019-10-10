import os
import shutil
import matlab.engine
import argparse
import datetime
import pandas as pd
import subprocess
import sys

from multiprocessing import Process
from multiprocessing import Pool

def ParsePQR(fileName):
    try:
        num_atom = 0;
        with open(fileName, "r") as f:

            tmp_arr = [];

            for l in f:
                if("ATOM" in l):
                    l = l.replace("\n", "");
                    l = l.replace("\t", " ");
                    tmp_arr.append(l.split());
                    num_atom += 1;

        col_names = ["Type", "Serial", "Name", "ChainID", "ResNum",
            "X", "Y", "Z", "Charge", "Radii"];

        df = pd.DataFrame(tmp_arr, columns=col_names);

        return df, num_atom;

    except FileNotFoundError as e:
        print(e);
    except InterruptedError as e:
        print(e);

def run_replica(apbsid, Nrandomdim, dfpqr):

    dirName = os.path.join("PyDsrar/output/Replicas","replica_"+str(apbsid));

    if(os.path.isdir(dirName)):
        shutil.rmtree(dirName);

    os.mkdir(dirName);

    #dfpqr, _ = ParsePQR(pqrName);

    xbase = "PyDsrar/output/APBS/x_sample_";
    ybase = "PyDsrar/output/APBS/y_sample_";
    zbase = "PyDsrar/output/APBS/z_sample_";

    dfx = pd.read_csv(xbase+str(apbsid)+"_list.dat", sep=" ", header=None);
    dfy = pd.read_csv(ybase+str(apbsid)+"_list.dat", sep=" ", header=None);
    dfz = pd.read_csv(zbase+str(apbsid)+"_list.dat", sep=" ", header=None);

    testmol_head = dfpqr[["Type", "Serial", "Name", "ChainID", "ResNum"]];

    qr_list  = dfpqr[["Charge", "Radii"]];

    locenergy = dirName+"/local_energy_" + str(apbsid)+".txt"
    gloenergy = dirName+"/global_energy_" + str(apbsid)+".txt"
    #os.path.touch(locenergy);
    #os.path.touch(gloenergy);

    for j in range(0, dfx.shape[1]):
        tx = dfx[j];
        ty = dfy[j];
        tz = dfz[j];

        tmp = pd.concat([testmol_head,tx,ty,tz,qr_list], axis=1);
        tmp.to_csv(dirName+"/testmol.pqr", sep=" ", header=False, index=False);

        with open(dirName+"/input", "w") as iw:
            with open("PyDsrar/input", "r") as ir:
                for l in ir:
                    if("FOLDERNAME" in l):
                        l = l.replace("FOLDERNAME", dirName);
                    iw.write(l);

        sys_str = "apbs "+dirName+"/input >" + dirName+ "/apbs_log.txt";
        #sys_str = "srun -N 1 -t 30:00 -p short -A STOCHASTIC_MOL ./apbs "+dirName+"/input >" + dirName+ "/apbs_log.txt";

        subprocess.call(sys_str, shell=True);

        with open(dirName+"/apbs_log.txt", "r") as r:
            for l in r:
                if("Local" in l):
                    strarr = l.split();
                    with open(locenergy, "a") as f:
                        f.write(strarr[-2]);
                        f.write("\n");
                elif("Global" in l):
                    strarr = l.split();
                    with open(gloenergy, "a") as f:
                        f.write(strarr[-2]);
                        f.write("\n");


if(__name__ == "__main__"):

    # parse command line arguments
    args = argparse.ArgumentParser()
    args.add_argument("--fileName", required=True, type=str, help="Molecule text file.");
    args.add_argument("--pqrName", required=True, type=str, help="PQR file name.");
    args.add_argument("--polyOrder", type=int, default=3, help="Required to construct orthogonal basis. Default is 3.");
    args.add_argument("--Nrandomdim", type=int, default=20, help="Number of random dimensions.");
    args.add_argument("--Npart", type=int, default=1000, help="Number of partitions. Default is 1000.");
    args.add_argument("--Nperpart", type=int, default=125, help="Number of columns per partition.");
    args.add_argument("--startN", type=int, default=200, help="Starting point for sampling in compute surrogate. Default is 200.");
    args.add_argument("--stopN", type=int, default=1400, help="End point for sampling in compute surrogate. Default is 1400. Cannot be larger than Npart * Nperpart");
    args.add_argument("--stepS", type=int, default=200, help="Step size for sampling in compute surrogate. Default is 200.");
    args.add_argument("--procs", type=int, default=2, help="Number of processors to use during the APBS caluclations. It may consume large amounts of memory.");

    p = args.parse_args();

    if(not os.path.isdir("PyDsrar/output/APBS")):
        os.mkdir("PyDsrar/output/");
        os.mkdir("PyDsrar/output/APBS");

    # init matlab engine
    eng = matlab.engine.start_matlab();
    eng.addpath("PyDsrar/");
    eng.addpath("PyDsrar/spgl1-1.9");

    Nrandomdim = p.Nrandomdim;
    Npart = p.Npart;
    polyOrder = p.polyOrder;
    Nperpart = p.Nperpart;

    dfpqr, Natom = ParsePQR(p.pqrName);


    # Step 1: Takes runs construct_MC_sample.m with input given in the command line.
    print("\n\nReading file: %s" % p.fileName);
    eng.construct_MC_sample(p.fileName, Nrandomdim, Npart, Nperpart, Natom, nargout=0);

    # Step 2: Perform APBS calculations
    dt0 = datetime.datetime.now();
    log_file_name = "log_"+str(dt0.month)+"_"+str(dt0.day)+"_";
    log_file_name += "_"+str(dt0.hour)+"_"+str(dt0.minute)+"_"+str(dt0.second);

    with open(os.path.join("PyDsrar/output", log_file_name), "w+") as logf:
        logf.write("Running time");
        res1 = datetime.datetime.now().timestamp();

        print("-- Running APBS");
        sys.stdout.flush();

        replicas_dir = os.path.join("PyDsrar", "output/Replicas");
        if(not os.path.isdir(replicas_dir)):
            os.mkdir(replicas_dir);

        ins = [];
        # TODO: change this line from 3 to the number of partitions (Npart)
        for i in range(1,Npart+1):
            ins.append((i, Nrandomdim, dfpqr));

        pool = Pool(processes=p.procs);
        pool.starmap(run_replica, ins);
        pool.close();
        pool.join();

        dt1 = datetime.datetime.now();
        dt = dt1 - dt0;
        dd  = dt.days;
        dh  = dt.seconds // 3600;
        dm  = (dt.seconds // 60)  % 60;

        tmpstr = "Total runtime: %d:%2d:%2d" % (dd, dh, dm);
        print(tmpstr);
        logf.write(tmpstr);
        logf.write("Finished at %s" % dt1.strftime("%m/%d/%y,%H:%M:%S"));

        mergeGlobal = "PyDsrar/output/global_energy_together.txt";
        if(os.path.isfile(mergeGlobal)):
            os.remove(mergeGlobal);

        with open(mergeGlobal, "a+") as fw:
            for i in range(1, Npart+1):
                strname = "PyDsrar/output/Replicas/";
                strname += "replica_" + str(i) + "/";
                strname += "global_energy_";
                strname += str(i) + ".txt";

                with open(strname, "r") as fr:
                    shutil.copyfileobj(fr, fw);

    # Step 3a
    eng.orthogonal_basis_construct(polyOrder, Nrandomdim, "random_MD_Traj.dat", nargout=0);

    # Step 3b
    eng.coeff_analysis_single(polyOrder, "coeff_scale_1.dat", nargout=0);

    # Step 4
    eng.compute_surrogate(polyOrder, Nrandomdim, "random_MD_Traj.dat",
    "global_energy_together.txt", p.startN, p.stopN, p.stepS,
    nargout=0);

    # Step 5, 6
    for i in range(p.startN, p.stopN+1, p.stepS):
        print("Running gradient_matrix_evaluation w/"+str(i));
        eng.gradient_matrix_evaluation_fun(i ,polyOrder, Nrandomdim, "random_MD_Traj.dat", nargout=0);

        print("-- orthogonal_basis_construct");
        eng.orthogonal_basis_construct(polyOrder, Nrandomdim, 'rotate_orth_basis_by_train_sam' + str(i) + '_rand_sample.dat', nargout=0);

        print("-- coeff_analysis_single")
        eng.coeff_analysis_single(polyOrder, "coeff_scale_1.dat", nargout=0);

        print('-- compute_surrogate');
        eng.compute_surrogate(polyOrder, Nrandomdim, "rotate_orth_basis_by_train_sam" + str(i) + "_rand_sample.dat",
        "global_energy_together.txt",i,i,p.stepS,
        nargout=0);
