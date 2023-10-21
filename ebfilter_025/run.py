#! /usr/bin/env python

import get_eb_score, control_count
import sys, os, re,  random
import pysam
from multiprocessing import Pool
from contextlib import closing
import process_samtools

def getvarcount_pysam(path_bam, targetMutationFile, n_mutation, mapping_qual_thres, base_qual_thres, filter_flags, verbose):

    hIN = open(targetMutationFile, 'r')
    
    h_bamfile = pysam.Samfile(path_bam, "rb")
        
    list_vardepth = [None] * n_mutation

    for k, line in enumerate(hIN):
        F = line.rstrip('\n').split('\t')
        chr, pos, _, ref, alt = F[0], int(F[1]), F[2], F[3], F[4]
        
        vardepths = control_count.varcount(
                        chr, pos, ref, alt, h_bamfile,
                        mapping_qual_thres, base_qual_thres, filter_flags, verbose)
                        
        list_vardepth[k] = vardepths

    hIN.close()
    h_bamfile.close()
    
    return list_vardepth
    
def getvarcount_samtools(pileup, targetMutationFile, n_mutation, mapping_qual_thres, base_qual_thres, verbose):

    list_mutation = [None] * n_mutation
    list_vardepth = [[0]*14] * n_mutation
    
    with open(targetMutationFile, 'r') as hIN:
        for k, line in enumerate(hIN):
            F = line.rstrip('\n').split('\t')
            chr, pos, _, ref, alt = F[0], int(F[1]), F[2], F[3], F[4]
            if alt == "-": pos -= 1
            
            list_mutation[k] = (chr, pos, ref, alt)

    with open(pileup) as f:
        for line in f:
            F = line.rstrip('\n').split('\t')
            chr, pos = F[0], int(F[1])
            
            k_muts = [k for k, x in enumerate(list_mutation) if chr == x[0] and pos == x[1]]
            
            for k_mut in k_muts:
                ref, alt = list_mutation[k_mut][2:]
                depth, baseBar, qual_list, mapqs, qnames = F[3:8]
                
                vardepths = control_count.varCount_samtools(
                                        ref, alt, 
                                        depth, baseBar, qual_list, mapqs, qnames,
                                        base_qual_thres, mapping_qual_thres, verbose)
                
                list_vardepth[k_mut] = vardepths
    
    return list_vardepth
    
def poolfun(args):
    fun = args[0]
    res = fun(*args[1:])
    return res
        
def EBFilter_worker(    targetMutationFile, 
                        targetBamPath, 
                        controlBamPathList, 
                        outputPath, 
                        mapping_qual_thres, 
                        base_qual_thres, 
                        filter_flags, 
                        region, 
                        method,
                        thread_num, 
                        epsilon, 
                        debug_mode):

    with open(controlBamPathList, 'r') as hcontrolbam:
        bams_ref = [x.strip('\n') for x in hcontrolbam.readlines() if x != '\n']
    
    n_mutation = sum([1 for _ in open(targetMutationFile)])
    n_reference = len(bams_ref)
    
    if method == 'pysam':
        # Randomizing the order to distribute access to the bam files was pointless.
        order_shuffleref = range(n_reference)
        # random.shuffle(order_shuffleref)
        
        list_vardepth_target = getvarcount_pysam(targetBamPath, targetMutationFile, n_mutation,
                                           mapping_qual_thres, base_qual_thres, filter_flags, False)
        
        list_list_vardepth_controls = [None] * n_reference
        
        if thread_num == 1:
            for k in order_shuffleref:
                bam_ref  = bams_ref[k]
                list_list_vardepth_controls[k] = getvarcount_pysam(bam_ref, targetMutationFile, n_mutation, 
                                                             mapping_qual_thres, base_qual_thres, filter_flags, True)
        else:
            args = [(getvarcount_pysam,
                    bams_ref[k], targetMutationFile, n_mutation,
                    mapping_qual_thres, base_qual_thres, filter_flags, True)
                    for k in order_shuffleref]
                    
            with closing(Pool( processes = thread_num-1 )) as p:
                    ress = p.map(poolfun, args)
                    p.terminate()
                    
            for res, k in zip(ress, order_shuffleref):
                list_list_vardepth_controls[k] = res

    else: #method == 'samtools':        
        ##########
        # generate pileup files
        process_samtools.anno2pileup(targetMutationFile, outputPath + '.pileup', targetBamPath, controlBamPathList, filter_flags, region, thread_num)
        ##########

        list_vardepth_target = getvarcount_samtools(
                outputPath + '.pileup_target', targetMutationFile, n_mutation, mapping_qual_thres, base_qual_thres, False)
        
        list_list_vardepth_controls = [None] * n_reference
        order_shuffleref = range(n_reference)

        if thread_num == 1:
            for k in order_shuffleref:
                pileup_file = outputPath + '.pileup_' + str(k)
                list_list_vardepth_controls[k] = getvarcount_samtools(
                        pileup_file, targetMutationFile, n_mutation, mapping_qual_thres, base_qual_thres, True)
        else:
            args = [(getvarcount_samtools,
                    outputPath + '.pileup_' + str(k), targetMutationFile, n_mutation, mapping_qual_thres, base_qual_thres, True)
                    for k in order_shuffleref]
                    
            with closing(Pool( processes = thread_num-1 )) as p:
                    ress = p.map(poolfun, args)
                    p.terminate()
                    
            for res, k in zip(ress, order_shuffleref):
                list_list_vardepth_controls[k] = res
                

    hIN = open(targetMutationFile, 'r')    
    hOUT = open(outputPath, 'w')
    
    for k, line in enumerate(hIN):
        F = line.rstrip('\n').split('\t')
        
        # If depths are 0 in all references, output is 0
        EB_score = 0 # if the variant is complex, we ignore that
        EB_score_whole = 0
        maxvarNum_ctrl  = "---"
        maxvarRate_ctrl = "---"
        maxvarNum_ctrl_whole   = "---"
        maxvarRate_ctrl_whole  = "---"
        maxvarNum_ctrl_strict  = "---"
        maxvarRate_ctrl_strict = "---"
        
        vardepth_target   = list_vardepth_target[k]
        vardepth_controls = [list_vardepth_control[k] for list_vardepth_control in list_list_vardepth_controls]
        
        if any([sum(x)>0 for x in vardepth_controls]):
            EB_score       = get_eb_score.get_eb_score(vardepth_target, vardepth_controls, epsilon)
            EB_score_whole = get_eb_score.get_eb_score(vardepth_target, vardepth_controls, epsilon, whole=True)
        
            varNum_ctrl   = [x[4] for x in vardepth_controls]
            depthNum_ctrl = [x[5] for x in vardepth_controls]
            varNum_ctrl_whole   = [x[10] for x in vardepth_controls]
            depthNum_ctrl_whole = [x[11] for x in vardepth_controls]
            varNum_ctrl_strict  = [x[12] for x in vardepth_controls]
            # varNum_ctrl_strict_whole = [x[13] for x in vardepth_controls]
        
            maxvarNum_ctrl  = max(varNum_ctrl)
            maxvarRate_ctrl = round( max([float(v) / d if d != 0 else 0 for v, d in zip(varNum_ctrl, depthNum_ctrl)]), 3)
            maxvarNum_ctrl_whole  = max(varNum_ctrl_whole)
            maxvarRate_ctrl_whole = round( max([float(v) / d if d != 0 else 0 for v, d in zip(varNum_ctrl_whole, depthNum_ctrl_whole)]), 3)
            maxvarNum_ctrl_strict  = max(varNum_ctrl_strict)
            maxvarRate_ctrl_strict = round( max([float(v) / d if d != 0 else 0 for v, d in zip(varNum_ctrl_strict, depthNum_ctrl)]), 3)
            
        # add the score and write the vcf record
        if debug_mode == False:    
            print >> hOUT, '\t'.join(F + [str(EB_score),       str(maxvarNum_ctrl),       str(maxvarRate_ctrl),
                                          str(EB_score_whole), str(maxvarNum_ctrl_whole), str(maxvarRate_ctrl_whole),
                                          str(maxvarNum_ctrl_strict), str(maxvarRate_ctrl_strict)
                                          ])
        else:
            print >> hOUT, '\t'.join(F + [str(EB_score), str(maxvarNum_ctrl), str(maxvarRate_ctrl),
                                          str(EB_score_whole), str(maxvarNum_ctrl_whole), str(maxvarRate_ctrl_whole),
                                          str(maxvarNum_ctrl_strict), str(maxvarRate_ctrl_strict),
                                          '{}/{}'.format(vardepth_target[4],vardepth_target[5]),
                                          ','.join(['{}/{}'.format(x,y) for x, y in zip(varNum_ctrl, depthNum_ctrl)])
                                     ])

    hIN.close()
    hOUT.close()
    

    # delete intermediate files
    if debug_mode == False and method == 'samtools':    
        os.remove(outputPath + '.pileup_target')
        
        for k in range(n_reference):
            os.remove(outputPath + '.pileup_' + str(k))


def main(args):

    # should add validity check for arguments
    targetMutationFile  = args.targetMutationFile
    targetBamPath       = args.targetBamPath
    controlBamPathList  = args.controlBamPathList
    outputPath          = args.outputPath

    mapping_qual_thres  = args.q
    base_qual_thres     = args.Q
    filter_flags        = args.ff
    thread_num          = args.t
    is_anno             = True if args.f == 'anno' else False
    is_loption          = args.loption
    region              = args.region
    epsilon             = args.epsilon
    method              = args.method
    debug_mode          = args.debug
    
    # file existence check
    if not os.path.exists(targetMutationFile):
        print >> sys.stderr, "No target mutation file: " + targetMutationFile
        sys.exit(1)

    if not os.path.exists(targetBamPath):
        print >> sys.stderr, "No target bam file: " + targetBamPath
        sys.exit(1)

    if not os.path.exists(targetBamPath + ".bai") and not os.path.exists(re.sub(r'bam$', "bai", targetBamPath)):
        print >> sys.stderr, "No index for target bam file: " + targetBamPath
        sys.exit(1)


    if not os.path.exists(controlBamPathList):
        print >> sys.stderr, "No control list file: " + controlBamPathList 
        sys.exit(1)

    with open(controlBamPathList) as hIN:
        for file in hIN:
            file = file.rstrip()
            if not os.path.exists(file):
                print >> sys.stderr, "No control bam file: " + file 
                sys.exit(1)

            if not os.path.exists(file + ".bai") and not os.path.exists(re.sub(r'bam$', "bai", file)):
                print >> sys.stderr, "No index control bam file: " + file 
                sys.exit(1)

    # single: samtools
    # pair  : pysam
    # Discriminate by 'bases / bases_tumor' column.
    if not is_loption or not region:
        method = 'pysam'
        
    elif method not in ('samtools', 'pysam'):
        with open(targetMutationFile, 'r') as f:
            firstitems = f.readline().strip('\n').split('\t')
            
            n_col_bases = 100
            for k, item in enumerate(firstitems):
                if ',' in item:
                    n_col_bases = k
                    break
                    
            if n_col_bases <= 8 : #single: 7
                method = 'samtools'
            else:
                method = 'pysam'

    if is_anno == True:
        EBFilter_worker(    targetMutationFile, 
                            targetBamPath, 
                            controlBamPathList, 
                            outputPath, 
                            mapping_qual_thres, 
                            base_qual_thres, 
                            filter_flags, 
                            region, 
                            method,
                            thread_num, 
                            epsilon, 
                            debug_mode)
        
    else: 
        print >> sys.stderr, "This version of EBfilter does not support vcf file."
        sys.exit(1)