#! /usr/bin/env python

import re, sys

ReIndel = re.compile('([\+\-])([0-9]+)([ACGTNacgtn]+)')
ReStart = re.compile('\^.')
ReEnd = re.compile('\$')

# verbose=True for control samples,
# in which indel is counted even when the length of indel is not exactly matched.
def varcount(chr, pos, ref, alt, h_bamfile, mapping_qual_thres, base_qual_thres, filter_flags, verbose):
#pos: 1-indexed = genomon output format

    exclude_flags = 0

    filter_flags_list = filter_flags.split(',')
    if 'UNMAP' in filter_flags_list:
        exclude_flags += 4
    if 'MUNMAP' in filter_flags_list:
        exclude_flags += 8
    if 'SECONDARY' in filter_flags_list:
        exclude_flags += 256
    if 'QCFAIL' in filter_flags_list:
        exclude_flags += 512
    if 'DUP' in filter_flags_list:
        exclude_flags += 1024
    if 'SUPPLEMENTARY' in filter_flags_list:
        exclude_flags += 2048
            
    #deletion
    if alt == '-':
        pos_0ind = pos - 2
    else:
        pos_0ind = pos - 1
    
    #flag_filter cannot remove SUPPLEMENTARY
    pileups = [x.pileups for x in h_bamfile.pileup(
                                    chr, pos_0ind, pos_0ind+1,
                                    truncate=True,
                                    ignore_overlaps=False,
                                    max_depth=100000,
                                    flag_filter=exclude_flags,
                                    flag_require=2,
                                    min_base_quality=0,
                                    min_mapping_quality=0)]
                                    
    if len(pileups) == 0:
        return [0] * 14
    else:
        pileups = pileups[0]
    
    # Reads with deletion over the base are excluded, because x.query_position is None in deletion.
    
    #SNV
    if ref != '-' and alt != '-':
        queries = [{'base':x.alignment.query_sequence[x.query_position],
                    'strand': '+' if not x.alignment.is_reverse else '-',
                    'qname':x.alignment.qname,
                    'mapq': x.alignment.mapping_quality,
                    'baseq': x.alignment.query_qualities[x.query_position]}
                    for x in pileups
                    if x.query_position is not None and 
                       (x.alignment.flag & 2) > 0 and 
                       (x.alignment.flag & exclude_flags) == 0]
                       
        queries_overthres = \
            [x for x in queries if x['mapq'] >= mapping_qual_thres and x['baseq'] >= base_qual_thres]
                       
        varreads_p_whole = [x for x in queries if (x['strand']=='+') and (x['base'] == alt) ]
        varreads_m_whole = [x for x in queries if (x['strand']=='-') and (x['base'] == alt) ]
        
        varreads_p = [x for x in varreads_p_whole if x in queries_overthres]
        varreads_m = [x for x in varreads_m_whole if x in queries_overthres]
        # varreads_all =len(set([ x['qname'] for x in queries_overthres if x['base'] == alt])) 
        if verbose:
            varreads_verbose_p = varreads_p
            varreads_verbose_m = varreads_m
            varreads_verbose_p_whole = varreads_p_whole
            varreads_verbose_m_whole = varreads_m_whole
            
    #indel
    else:
                    
        #insertion
        if ref == '-':
            queries = [{'indel':x.indel,
                        'strand': '+' if not x.alignment.is_reverse else '-',
                        'qname':x.alignment.qname,
                        'seq':x.alignment.query_sequence,
                        'query_pos':x.query_position,
                        'is_same': False,
                        'mapq': x.alignment.mapping_quality}
                        for x in pileups
                        if x.query_position is not None and
                           (x.alignment.flag & 2) > 0 and
                           (x.alignment.flag & exclude_flags) == 0]# and 
                            #x.alignment.query_qualities[x.query_position] >= base_qual_thres]
                            
                            #quality is incorporated for compatibility
                        
            if verbose:
                varreads_verbose_p_whole = \
                    [x for x in queries if (x['strand']=='+') and (x['indel'] > 0) ]
                varreads_verbose_m_whole = \
                    [x for x in queries if (x['strand']=='-') and (x['indel'] > 0) ]
                # varreads_all =len(set([ x['qname'] for x in queries if x['indel'] > 0])) 

            for k, query in enumerate(queries):
                if query['indel'] == len(alt):
                    ins_bases = query['seq'][ query['query_pos']+1 : query['query_pos']+query['indel']+1]
                    if ins_bases == alt:
                        queries[k]['is_same'] = True
        
            varreads_p_whole = [x for x in queries if (x['strand']=='+') and x['is_same'] ]
            varreads_m_whole = [x for x in queries if (x['strand']=='-') and x['is_same'] ]
                
                # varreads_all =len(set([ x['qname'] for x in queries if x['is_same']])) 
        
        #deletion
        else:
            queries = [{'indel':x.indel,
                        'strand': '+' if not x.alignment.is_reverse else '-',
                        'qname':x.alignment.qname,
                        'mapq': x.alignment.mapping_quality}
                        for x in pileups
                        if x.query_position is not None and
                           (x.alignment.flag & 2) > 0 and
                           (x.alignment.flag & exclude_flags) == 0]#and 
                            #x.alignment.query_qualities[x.query_position] >= base_qual_thres]
                            
                            #quality is incorporated for compatibility
                            
            if verbose:
                varreads_verbose_p_whole = \
                    [x for x in queries if (x['strand']=='+') and (x['indel'] < 0) ]
                varreads_verbose_m_whole = \
                    [x for x in queries if (x['strand']=='-') and (x['indel'] < 0) ]
                # varreads_all =len(set([ x['qname'] for x in queries if x['indel'] < 0])) 

            varreads_p_whole = [x for x in queries if (x['strand']=='+') and (x['indel'] == -len(ref)) ]
            varreads_m_whole = [x for x in queries if (x['strand']=='-') and (x['indel'] == -len(ref)) ]
            # varreads_all =len(set([ x['qname'] for x in queries if x['indel'] == -len(ref)])) 

        #indel
        queries_overthres = [x for x in queries if x['mapq'] >= mapping_qual_thres]
        varreads_p = [x for x in varreads_p_whole if x in queries_overthres]
        varreads_m = [x for x in varreads_m_whole if x in queries_overthres]
        if verbose:
            varreads_verbose_p = [x for x in varreads_verbose_p_whole if x in queries_overthres]
            varreads_verbose_m = [x for x in varreads_verbose_m_whole if x in queries_overthres]
            
    depth_p   = sum([ x['strand']=='+' for x in queries_overthres])
    depth_m   = sum([ x['strand']=='-' for x in queries_overthres])
    depth_all = len(set([ x['qname'] for x in queries_overthres])) 
    
    depth_p_whole = sum([ x['strand']=='+' for x in queries])
    depth_m_whole = sum([ x['strand']=='-' for x in queries])
    depth_all_whole = len(set([ x['qname'] for x in queries ]))
    
    if verbose:
        num_varread_p_whole = len(varreads_verbose_p_whole)
        num_varread_m_whole = len(varreads_verbose_m_whole)
        num_varread_p       = len(varreads_verbose_p)
        num_varread_m       = len(varreads_verbose_m)
        
        num_varread_all_whole = len(set([x['qname']
                    for x in varreads_verbose_p_whole + varreads_verbose_m_whole]))
        num_varread_all       = len(set([x['qname']
                    for x in varreads_verbose_p + varreads_verbose_m]))
        
        num_varread_strict_all_whole = len(set([x['qname']
                    for x in varreads_p_whole + varreads_m_whole]))
        num_varread_strict_all       = len(set([x['qname'] for x in varreads_p + varreads_m]))
            
        return (num_varread_p,   depth_p,
                num_varread_m,   depth_m,
                num_varread_all, depth_all,
                num_varread_p_whole,   depth_p_whole,
                num_varread_m_whole,   depth_m_whole,
                num_varread_all_whole, depth_all_whole,
                num_varread_strict_all,
                num_varread_strict_all_whole)
                
    else:
        num_varread_p_whole = len(varreads_p_whole)
        num_varread_m_whole = len(varreads_m_whole)
        num_varread_p       = len(varreads_p)
        num_varread_m       = len(varreads_m)
        
        num_varread_all_whole = len(set([x['qname']
                    for x in varreads_p_whole + varreads_m_whole]))
        num_varread_all       = len(set([x['qname'] for x in varreads_p + varreads_m]))
        
        return (num_varread_p,   depth_p,
                num_varread_m,   depth_m,
                num_varread_all, depth_all,
                num_varread_p_whole,   depth_p_whole,
                num_varread_m_whole,   depth_m_whole,
                num_varread_all_whole, depth_all_whole,
                0, 0)
        
    
    # this is compatible with old, but wrong!
    # depth_all = len(set([ x['qname'] for x in queries ] + 
                        # [ x.alignment.qname for x in pileups
                        # if x.is_del and 
                       # (x.alignment.flag & 2) > 0 and 
                       # (x.alignment.flag & exclude_flags) == 0]))
            
        
        

# verbose=True for control samples,
# in which indel is counted even when the length of indel is not exactly matched.
def varCount_samtools(ref, alt, 
                      depth, baseBar, qual_list, mapqs, qnames,
                      base_qual_thres, mapping_qual_thres, verbose):
        
    var = ""
    if ref != "-" and alt != "-":
        var = alt
    else:
        if ref == "-":
            var = "+" + alt
        elif alt == "-":
            var = "-" + ref
            
    if depth == 0: return [0] * 14

    results = [0] * 14
    
    baseBar = ReStart.sub( '', baseBar )
    baseBar = baseBar.replace('$','')
    
    isindel = (ref == '-') | (alt == '-')
    
    qname_list = qnames.split(',')
    qnames_alt       = {'p':[], 'm':[]}
    qnames_alt_whole = {'p':[], 'm':[]}
    qnames_verbose_alt       = {'p':[], 'm':[]}
    qnames_verbose_alt_whole = {'p':[], 'm':[]}

    #
    # Look for deletion/insertion and save info in 'indel' dictionary
    #
    #   ([\+\-])[0-9]+[ACGTNacgtn]+
    #
    # m.group(1): + or - (deletion/insertion)
    # m.group(2): number of deletion/insertion
    # m.group(3): nucleotides
    #
    deleted = 0
    iter = ReIndel.finditer( baseBar )
    for m in iter:
        site = m.start()
        indeltype = m.group( 1 )
        num = m.group( 2 )
        varChar = m.group( 3 )[ 0:int( num ) ]
        
        k_base_before = site - deleted - 1
        baseqname = qname_list[ k_base_before ]
        isbasemapq = 33 + mapping_qual_thres <= ord(mapqs[ k_base_before ])

        if isindel and var[0] == indeltype:
            if indeltype == '+':
                is_same = var[1:] == varChar.upper()
            else:
                is_same = len(var[1:]) == len(varChar)
                
            isplus = 'p' if varChar.isupper() else 'm'
            
            if verbose:
                qnames_verbose_alt_whole[isplus].append(baseqname)
                if isbasemapq:
                    qnames_verbose_alt[isplus].append(baseqname)
                    
            if is_same:
                qnames_alt_whole[isplus].append(baseqname)
                if isbasemapq:
                    qnames_alt[isplus].append(baseqname)
                            
        baseBar = baseBar[ 0:site - deleted ] + baseBar[ site + int(num) + len( num ) + 1 - deleted: ]
        deleted += 1 + len( num ) + int( num )
        
    # error check
    if any([ len(baseBar) != len(qual_list),
             len( baseBar ) != len(mapqs),
             len( baseBar ) != len(qname_list) ]):
        print >> sys.stderr, baseBar + '\n' + qual_list + mapqs
        print >> sys.stderr, "lengths of bases and qualities are different!"
        sys.exit(1)
        
    #
    # Count mismatch
    #

    isplus  = [x.isupper() for x in baseBar]
    isminus = [x.islower() for x in baseBar]
    isnotdel = [(x != '#') & (x != '*') for x in baseBar]
    # exclude # and *
    
    ismapqual   = [33 + mapping_qual_thres <= ord(mapq) for mapq in mapqs]
    isqual      = [(33 + base_qual_thres <= ord(baseq)) & (33 + mapping_qual_thres <= ord(mapq))
                    for baseq, mapq in zip(qual_list, mapqs)]
    
    if isindel:
        depth_p = sum([x & y for x, y in zip(ismapqual,  isplus  )])
        depth_m = sum([x & y for x, y in zip(ismapqual,  isminus )])
        depth_all = len(set([qname for qname, x,y in zip(qname_list, isnotdel, ismapqual) if x & y]))
        if verbose:
            results[0] = len(qnames_verbose_alt['p'])
            results[2] = len(qnames_verbose_alt['m'])
            results[4] = len(set(qnames_verbose_alt['p'] + qnames_verbose_alt['m']))
            results[6] = len(qnames_verbose_alt_whole['p'])
            results[8] = len(qnames_verbose_alt_whole['m'])
            results[10]= len(set(qnames_verbose_alt_whole['p'] + qnames_verbose_alt_whole['m']))
            results[12]= len(set(qnames_alt['p'] + qnames_alt['m']))
            results[13]= len(set(qnames_alt_whole['p'] + qnames_alt_whole['m']))
        else:
            results[0] = len(qnames_alt['p'])
            results[2] = len(qnames_alt['m'])
            results[4] = len(set(qnames_alt['p'] + qnames_alt['m']))
            results[6] = len(qnames_alt_whole['p'])
            results[8] = len(qnames_alt_whole['m'])
            results[10]= len(set(qnames_alt_whole['p'] + qnames_alt_whole['m']))
    else:
        depth_p   = sum([x & y for x, y in zip(isqual,  isplus )])
        depth_m   = sum([x & y for x, y in zip(isqual,  isminus)])
        depth_all = len(set([qname for qname, q1, q2 in zip(qname_list, isqual, isnotdel) if q1 and q2]))
        results[0] = sum([q & (base == var)         for base, q in zip(baseBar, isqual)])
        results[2] = sum([q & (base == var.lower()) for base, q in zip(baseBar, isqual)])
        results[4] = len(set([qname for base, q, qname in zip(baseBar, isqual, qname_list) if q and (base.upper() == var)]))
        results[6] = baseBar.count(var)
        results[8] = baseBar.count(var.lower())
        results[10]= len(set([qname for base, qname in zip(baseBar, qname_list) if base.upper() == var]))
        results[12]= results[4]
        results[13]= results[10]
        
    depth_p_whole   = sum(isplus)
    depth_m_whole   = sum(isminus)
    depth_all_whole = len(set([qname for qname, x in zip(qname_list, isnotdel) if x]))
    
    results[1] = depth_p
    results[3] = depth_m
    results[5] = depth_all
    results[7] = depth_p_whole
    results[9] = depth_m_whole
    results[11] = depth_all_whole
    
    return results