#! /usr/bin/env python

import sys, os, subprocess, time
        
def anno2pileup(inputFilePath, outputFilePath, targetBamPath, bamPath, filter_flags, region, thread_num):

    hBAM = open(bamPath, 'r')
    bams_ref = [x.strip('\n') for x in hBAM.readlines() if x != '\n']
    bams = [targetBamPath] + bams_ref
    suffixes = ['_target'] + ['_'+str(k) for k in range(len(bams_ref))]
    hBAM.close()


    hIN = open(inputFilePath, 'r')
    
    # make bed file for mpileup
    hOUT2 = open(outputFilePath + ".region_list.bed", 'w')
    for k, line in enumerate(hIN):
    
        F = line.rstrip('\n').split('\t')
        if F[4] == "-": # for deletion in anno format
            print >> hOUT2, F[0] + '\t' + str(int(F[1]) - 2) + '\t' + str(int(F[1]) - 1) 
        else:
            print >> hOUT2, F[0] + '\t' + str(int(F[1]) - 1) + '\t' + F[1]

    hOUT2.close()

    ### Fix:0.2.2 add '-x' to exactly count bases for overlap read.
    samtools_mpileup_command = \
        ["samtools", "mpileup", "--output-QNAME", "-B", "-x", "-d", "10000000",
         "-q", "0", "-Q", "0", "-s", "--ff", filter_flags,
         "-l", outputFilePath + ".region_list.bed"]

    if region != "":
        samtools_mpileup_command = samtools_mpileup_command + ["-r", region]
        
    FNULL = open(os.devnull, 'w')
    hOUTs = []
    procs = []
    for suffix, bam in zip(suffixes, bams):
    
        hOUT = open(outputFilePath + suffix, 'w')
        hOUTs.append(hOUT)
        
        samtools_mpileup_command_bam = samtools_mpileup_command + [bam]

        if thread_num <= 1:
            proc = subprocess.check_call(samtools_mpileup_command_bam, 
                                         stdout = hOUT, stderr = FNULL)
        # multi-threading mode
        else:
            proc = subprocess.Popen(samtools_mpileup_command_bam, 
                                         stdout = hOUT, stderr = FNULL)
            procs.append(proc)
            
            while True:
                if sum([x.poll() is None for x in procs]) >= thread_num-1:
                    time.sleep(1)
                else:
                    break
    
    if thread_num >= 2:
        for proc in procs:
            proc.wait()

    for hOUT in hOUTs:
        hOUT.close()

    os.remove(outputFilePath + ".region_list.bed")
    hIN.close()
    FNULL.close()
    