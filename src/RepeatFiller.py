#!/usr/bin/env python

# Ekaterina Osipova, MPI-CBG/MPI-PKS, 2018

##########

import sys
import os
import subprocess
import argparse
import logging
import re
import io
import time
import tempfile
############################


############################
# Makes a list of jobs to run in temp shell script
#
# inChainFile: file containing all chains
# outf: path to output file; shell commands will be written into this file 
def MakeShellList(inChainFile, outf):

    # open output file or fail
    try:
        fhout = open(outf, 'w')
    except IOError as e:
        logging.error('Cannot write to shell script', e.errno, e.strerror)
        sys.exit(1)

    # write header for shell script
    fhout.write("#!/usr/bin/env bash\n")
    fhout.write("#set -o pipefail\n")
    fhout.write("#set -e\n")

    # Numbers for tracking the line number
    lineNmbr = 0

    # two bit files
    T2bit = args.T2bit
    Q2bit = args.Q2bit

    lastzVar = args.lastz
    axtChainVar = args.axtChain
    chainSortVar = args.chainSort

    # read inChainFile line by line
    with open(inChainFile, 'r') as inChain:

        try:
            for line in inChain:
                lineNmbr += 1

                ll = line.split()
                if len(ll) > 0:
                    if ll[0] == 'chain':
                        # read the chain header:
                        # e.g. chain 196228 chr4 62094675 + 12690854 12816143 chr23 24050845 - 20051667 20145391 1252
                        score = int(ll[1])
                        tName, tStart, tEnd = ll[2], int(ll[5]), int(ll[6])
                        qName, qStartx, qEndx, qStrand = ll[7], int(ll[10]), int(ll[11]), ll[9]
                        qSize = int(ll[8])
                        logging.info('qStrand = {}'.format(qStrand))
                        # changing coords for -strand if necessary

                        qStart = qStartx
                        qEnd = qEndx
                        lastzParameters = args.lastzParameters + ' --strand=plus'

                        if qStrand == "-":
                            lastzParameters = args.lastzParameters + ' --strand=minus'

                        if ll[4] != "+":
                            logging.error("ERROR: target strand is not + for chain:{}".format(line))
                            sys.exit(1)
                        # check if we consider this chain
                        logging.info("score of this chain = {}".format(score))
                        if ((score >= args.chainMinScore) and (tEnd - tStart >= args.chainMinSizeT) and (qEnd - qStart >= args.chainMinSizeQ)):
                            logging.info("valid chain")
                            curTPos = tStart
                            curQPos = qStart

                            line = next(inChain)
                            lineNmbr += 1
                            while(re.match("^\d+", line) is not None):
                                a = line.split()
                                if (len(a) == 1):
                                    TblockEnd = curTPos + int(a[0])
                                    QblockEnd = curQPos + int(a[0])
                                    logging.info("it was the last block\n")
                                else:
                                    blockLen = int(a[0])
                                    TblockEnd = curTPos + blockLen
                                    QblockEnd = curQPos + blockLen
                                    TgapEnd = curTPos + blockLen + int(a[1])
                                    QgapEnd = curQPos + blockLen + int(a[2])
                                    tGapSpan = TgapEnd - TblockEnd
                                    qGapSpan = QgapEnd - QblockEnd
                                    # check if we want to fill this gap
                                    if ((tGapSpan >= args.gapMinSizeT) and (tGapSpan <= args.gapMaxSizeT) and (
                                        qGapSpan >= args.gapMinSizeQ) and (qGapSpan <= args.gapMaxSizeQ)):
                                        logging.info("yes, this gap will be filled: {}".format(line.strip()))
                                        TblockEnd += 1
                                        QblockEnd += 1

                                        # replace the content of the unmask by '[unmask]' if the user DID NOT set the flag, otherwise ''
                                        if args.keepmask:
                                            unmask = ''
                                        else:
                                            unmask = '[unmask]'

                                        if qStrand == "-":
                                            realQblockEnd = qSize - QgapEnd + 1
                                            realQgapEnd = qSize - QblockEnd + 1
                                        else:
                                            realQblockEnd = QblockEnd
                                            realQgapEnd = QgapEnd

                                        logging.info("running lastz on the block:")
                                        regionToBeFilled = [tName, str(TblockEnd), str(TgapEnd), qName, str(realQblockEnd), str(realQgapEnd)]
                                        logging.info(' '.join(regionToBeFilled))

                                        # making lastz command for this region
                                        command1 = ' {0}/{1}[{2}..{3}]{9} {4}/{5}[{6}..{7}]{9} --format=axt {8} | '.format(T2bit, tName, TblockEnd, TgapEnd, Q2bit, qName, realQblockEnd, realQgapEnd, lastzParameters, unmask)
                                        command2 = ' -linearGap=loose stdin {0} {1} stdout 2> /dev/null | '.format(T2bit, Q2bit)
                                        command3 = ' stdin stdout'
                                        command_lastz = lastzVar + command1 + axtChainVar + command2 + chainSortVar + command3
                                        ### adding this lastz run to a shell command list; lineNmbr - 1 because we start with 1 and later with 0
                                        shellCommand = 'echo -e "LINE{0}\\n{1}\\n{2}\\n{3}\\n{4}\\n{5}\\n"; {6}; echo -e "LINE{0}\\n"\n'.\
                                            format((lineNmbr-1), blockLen, TblockEnd, TgapEnd, realQblockEnd, realQgapEnd, command_lastz)
                                        fhout.write(shellCommand)

                                    curQPos = QgapEnd
                                    curTPos = TgapEnd

                                # get next line, break if no more lines in string
                                try:
                                     line = next(inChain)
                                     lineNmbr += 1
                                except Exception:
                                     break
                        else:
                            logging.info("invalid chain\n")

                            # save chain header line; get next line
                            line = next(inChain)
                            lineNmbr += 1

                            # read the rest of the chain store blocks
                            while(re.match("^\d+", line) is not None):
                                try:
                                    line = next(inChain)
                                    lineNmbr += 1
                                except Exception:
                                    break
        except StopIteration:
            pass
    logging.info("Done with reading gaps")
    logging.info('\n')
    logging.info('\n')
    
    # close file handler
    fhout.close()

    return()

######################
# Takes temp file with all shell commands to run and returns lastz output in a single string
#
def RunAllShell(shellfile):
    AllShellCommand = 'bash ' + shellfile

    try:
        allMiniChains = subprocess.check_output(AllShellCommand, shell=True)

    except subprocess.CalledProcessError as shellrun:
        sys.exit(1)
    allMiniChains = allMiniChains.decode()

    return(allMiniChains)



######################
# Takes the whole lastz output chain list and returns list containing chain block starting at cur_position
# 
# (lineNmbr: [blockLen, TblockEnd, TgapEnd, realQblockEnd, realQgapEnd, 'allChainsStrings']) from LINE# to LINE#
# returns a dictionary
#
# allMiniChainsSplit: list containing line-wise lastz output file includung LINE### statements as block separator
# positions: starting position of current block (should start with LINE###
#
# returns list of line split strings: one output block starting with LINE### and ending with LINE#
# raises ValueError if block not properly separated by LINE#
#
def GetChainBlockFromLastzOutput(allMiniChainsSplit, cur_position):
    
    position = cur_position
    start = position
    end = None
    line = allMiniChainsSplit[position]
    reLine = re.compile(r'LINE\d+')

    # check whether initial line start with LINE#, raise error otherwise
    if reLine.match(line) is not None:
        
        position += 1 #get next line 
        line = allMiniChainsSplit[position]
        
        # process block until LINE# as end separator is encountered
        while reLine.match(line) is None:
              position += 1 #get next line 
              line = allMiniChainsSplit[position]

        # check that last line contains LINE#
        if reLine.match(line) is not None: 
            end = position
        else:
            raise ValueError('ERROR! allMiniChainsSplit end separator line at position ' + str(position) + ' does not start with LINE...')

    else:
        raise ValueError('ERROR! allMiniChainsSplit start separator line at position ' + str(position) + ' does not start with LINE...')

    curBlockList = allMiniChainsSplit[start:(end+1)]

    return(curBlockList)


######################
# Takes first chain from a chain list
# returns the header: "chain 52633 chr..." and a list of lines of this chain
# returns twice None if no chains are present
#
def TakeFirstChainFromList(chainList):
    
    headLine = None
    chainContent = None
    chainStart = None
    chainEnd = None

    for pos in range(0, len(chainList)):

        line = chainList[pos]

        # check if chain line
        m = re.match(r'chain', line)
        if m is not None:
            headLine = line.strip('\n')
            
            # process and store end position
            pos +=1
            line = chainList[pos]
            chainStart = pos
            
            while(re.match("^\d+", line) is not None):                
                pos +=1
                line = chainList[pos]            
            chainEnd = pos # actually position after chain            

            # don't process lower scoring chains
            break 

    # extract chain
    if chainStart is not None:
        chainContent = chainList[chainStart:chainEnd]

    return(headLine, chainContent)


###################
# After filling chain we need to insert it back on the right place
# InsertChainContent calculates new coordinates for a chain to be inserted
# and returns a list of lines, that were changed in comparision with an old chain file

def InsertChainContent(ChainContent, BestChain, BlockLenA, TblockEnd, TgapEnd, LoQblockEnd, LoQgapEnd):
    TlastzStart = int(BestChain.split()[5]) + 1
    TlastzEnd = int(BestChain.split()[6])

    if BestChain.split()[9] == '+':
        QlastzStart = int(BestChain.split()[10]) + 1
        QlastzEnd = int(BestChain.split()[11])
    else:
        # recalculate -strand coords to +strand:
        QlastzStart = int(BestChain.split()[8]) - int(BestChain.split()[10])
        QlastzEnd = int(BestChain.split()[8]) - int(BestChain.split()[11]) + 1

        tempQ = LoQgapEnd
        LoQgapEnd = LoQblockEnd
        LoQblockEnd = tempQ

    BlockToAdd = []

    if BestChain.split()[9] == '+':
        FirstQgap = abs(QlastzStart - int(LoQblockEnd))
        LastQGap = abs(int(LoQgapEnd) - QlastzEnd)
    else:
        FirstQgap = abs(QlastzStart - int(LoQblockEnd))
        LastQGap = abs(int(LoQgapEnd) - QlastzEnd)

    FirstLine = str(BlockLenA) + '\t' + str(TlastzStart - int(TblockEnd)) + '\t' + str(FirstQgap) + '\n'

    BlockToAdd.append(FirstLine)
    for i in range(0, len(ChainContent)-1):
        BlockToAdd.append(ChainContent[i])

    LastLine = ChainContent[len(ChainContent)-1].strip() + '\t' + str(int(TgapEnd)-TlastzEnd) +\
                   '\t' + str(LastQGap) + '\n'
    BlockToAdd.append(LastLine)
    return(BlockToAdd)



######## main ###########

# tracking runtime
start_time = time.time()


###### Optional arguments #######
parser = argparse.ArgumentParser(description='This script processes chain file, finds gaps and builds local alignments using lastz, then inserts new blocks to the original chain file', \
                                 epilog='Example of use: RepeatFiller.py -c hg38.speTri2.all.chain -T2 hg38.2bit -Q2 speTri2.2bit')
parser.add_argument('-v', '--verbose', action='store_true', help="if -v is not specified, only ERROR messages will be shown")
parser.add_argument('--lastz', '-l', type=str, default='lastz', help="path to lastz executable, default = lastz")
parser.add_argument('--axtChain', '-x', type=str, default='axtChain', help="path to axtChain executable, default = axtChain")
parser.add_argument('--axtChainParameters', '-xparam', type=str, default='linearGap=loose', help='parameters to run axtChain, default: linearGap=loose')
parser.add_argument('--chainSort', '-s', type=str, default='chainSort', help="path to chainSort executable, default = chainSort")
parser.add_argument('--output', '-o', type=str, help="name of output chain file. If not specified chains go to stdout")
parser.add_argument('--workdir', '-w', type=str, default='./', help="working directory for temp files, default = ./")

# initial parameters
parser.add_argument('--chainMinScore', '-mscore', type=int, default=25000, help="consider only chains with a chainMinScore, default mscore = 25000")
parser.add_argument('--chainMinSizeT', '-mst', type=int, default=0, help="consider only chains with a chainMinSizeT, default consider all")
parser.add_argument('--chainMinSizeQ', '-msq', type=int, default=0, help="consider only chains with a chainMinSizeQ, default consider all")
parser.add_argument('--gapMinSizeT', '-gmint', type=int, default=30, help="fill only gaps that are at least that long on the target side, default gmint = 30")
parser.add_argument('--gapMinSizeQ', '-gminq', type=int, default=30, help="fill only gaps that are at least that long on the query side, default gminq = 30")
parser.add_argument('--gapMaxSizeT', '-gmaxt', type=int, default=20000, help="fill only gaps that are at most that long on the target side, default gmaxt = 20000")
parser.add_argument('--gapMaxSizeQ', '-gmaxq', type=int, default=20000, help="fill only gaps that are at most that long on the query side, default gmaxq = 20000")
parser.add_argument('--lastzParameters', '-lparam', type=str, default=' K=2000 L=3000 M=0 T=0 W=5 ', help="line with lastz parameters, default 'K=2000 L=3000 M=0 T=0 W=5' ")
parser.add_argument('--keepmask', '-km', action='store_false', help="KEEPS masking (lower case characters) in 2bit files if specified")
parser.add_argument('--scoreThreshold', '-st', type=int, default=5000, help="insert only chains that have at least this score, default st = 5000")

###### Required arguments #######
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument("--chain", "-c", type=str, help="all.chain file", required=True)
requiredNamed.add_argument("--T2bit", '-T2', type=str, help="reference 2bit file", required=True)
requiredNamed.add_argument("--Q2bit", '-Q2', type=str, help="query 2bit file", required=True)
args = parser.parse_args()



if args.verbose:
   logging.basicConfig(level=logging.INFO)


# 1) Loop through .all.chain and make a jobList in a shell script
if not os.path.isdir(args.workdir):
    logging.error("ERROR! Working directory '" + args.workdir + "' does not exist.")
    sys.exit(1)

# create a file for writing jobList
try:
    temp = tempfile.NamedTemporaryFile(prefix = "tempCGFjobList", dir = args.workdir, delete = False )
    temp.close()
except PermissionError as e:
    logging.error('ERROR! Failed to create temporary file inside \'{}\'. {} {}'.format(args.workdir, e.errno, e.strerror) )
    sys.exit(1)

# Find gaps and write corresponding jobs to a shell script
MakeShellList(args.chain, temp.name)


# 2) Run the prepared shell script
allMiniChains = RunAllShell(temp.name)

# remove temp file
os.unlink(temp.name)

# string will be used for final output
OutputChain = ''

# position of next block in miniChain output 
next_pos = 0
# next line number where chain is filled
next_lineNmbr = None

# regexp for getting lineNumber
relineNmbr = re.compile(r'LINE(\d+)')

# list of mini Chain Blocks
allMiniChainLines = [i + '\n' for i in allMiniChains.split('\n')]
lenAllMiniChainLines = len(allMiniChainLines)

# get current mini chain
curMiniBlockLines = GetChainBlockFromLastzOutput(allMiniChainLines, next_pos)
# get next lineNmbr where intial chain will be filled with gaps
m = relineNmbr.match(curMiniBlockLines[0])
if m is not None:
    next_lineNmbr = int(m.group(1))
else:
    raise ValueError('ERROR! Could not extract line number from separator current miniChain block')



### process initial chain and fill gaps with miniChainBlocks

# write output; open file handle
if args.output:
    try:
        ouf = open(args.output, 'w')
    except IOError as e:
        logging.error('Cannot write to outputfile', e.errno, e.strerror)
        sys.exit(1)


with open(args.chain, 'r') as inChain:
    # initialize line counter in original chain file
    lineNmbr = 0

    try:
        for line in inChain:
            # empty output string
            OutputChain = ""
            # update chain
            if lineNmbr == next_lineNmbr:
                # strip first and last line containing LINE# from block
                valueList = curMiniBlockLines[1:(len(curMiniBlockLines) - 1)]
                Coords = valueList[:5]
                # remove new lines from Coords elements
                Coords = map(lambda s: s.strip(), Coords)

                # update next_pos and get next mini chain block; +1 since we have new line after each block in the output
                next_pos = next_pos + len(curMiniBlockLines) + 1
                # test that we are not out of bounds, i.e. last entry, -1 since last line is new line
                if next_pos < lenAllMiniChainLines - 1:
                    curMiniBlockLines = GetChainBlockFromLastzOutput(allMiniChainLines, next_pos)
                    # get next lineNmbr
                    m = relineNmbr.match(curMiniBlockLines[0])
                    if m is not None:
                        next_lineNmbr = int(m.group(1))
                    else:
                        raise ValueError('ERROR! Could not extract line number from separator current miniChain block')

                # get chain to be inserted
                BestChain, ChainContent = TakeFirstChainFromList(valueList[5:])

                # insert nothing if no chain in block
                if BestChain is not None:
                    if int(BestChain.split()[1]) >= args.scoreThreshold:
                        logging.info('Best lastz output chain = {}'.format(BestChain))

                        InsertBlock = InsertChainContent(ChainContent, BestChain, *Coords)
                        OutputChain = ''.join(InsertBlock)
                        logging.info('--- %s seconds ---' % (time.time() - start_time))
                    else:
                        logging.info("lastz output chains have low score\n")
                        OutputChain = line
                else:
                    logging.info("lastz changed nothing in this block\n")
                    OutputChain = line
            else:
                # Just add this line to the chain string
                OutputChain = line

            # go to the next line of the file
            lineNmbr += 1

            # print output
            if args.output:
                ouf.write(OutputChain)
            else:
                print(OutputChain)
    except StopIteration:
        pass

# close output file handle
if args.output:    
    ouf.close()

logging.info('--- Final runtime: %s seconds ---' % (time.time() - start_time))
