import re
import os
import numpy as np
from tools import reverse_seq

def load_position(file):

    f = open(file,'r')
    line = f.readline()

    fw_join = open(file+'_position.txt', 'w')

    while line:
        if line.startswith('>'):
            id = line.split(' ')[0]
            result = re.search(r'location=(.*)', line).group().replace('..>','..').replace('<','')
            if 'join' in result:
                position = result.split('join(')[1].split(')')[0]
                if 'complement' in result:
                    fw_join.write(id + '\t' + '1' + '\n')
                    fw_join.write(position + '\n')
                else:
                    fw_join.write(id + '\t' + '0' + '\n')
                    fw_join.write(position + '\n')
        line = f.readline()

    fw_join.close()

    f.close()


def load_seq(file):

    f = open(file, 'r')
    line = f.readline()

    fw_seq = open(file + '_seq.txt', 'w')

    while line:
        if line.startswith('>'):
            id = line.split(' ')[0]
            fw_seq.write(id+'\n')
        else:
            fw_seq.write(line)
        line = f.readline()


    fw_seq.close()

    f.close()


def load_position_info(positionfile):
    """
    Reads Alternetive Splicing position information
    :param positionfile: File path to position information
    :return: position dictionary
            key: AS id e.g. >lcl|NW_015379189.1_mrna_XM_015765349.1_50431
            value: [complementary marker, [int number position marker]] e.g. ['0', [3401, 4025, 4160, 4408, 4516, 5122, 5221, 5288, 5384, 5446, 5531, 5689]]
    """
    position = dict()

    f = open(positionfile, 'r')
    line = f.readline()
    while line:
        if line.startswith('>'):
            id = line.strip().split('\t')[0]
            comp = line.strip().split('\t')[1]
            tmp = [int(x) for x in f.readline().strip().replace(',', '..').split('..')]
            position[id] = [comp, tmp]
            # print(id)
            # print(position[id])
        line = f.readline()
    f.close()

    return position


def load_sequences_segment(seqdir):
    """
    Reads genomic segmented sequences
    :param seqdir: 8000bp genomic sequences segment file path
    :return: sequences dictionary
            key: genomic sequences segments id e.g. >NC_029260.1   1
            value: segment sequences e.g. 'cctaaaccctaaaccctaaaccctaaaccctaaacc'
    """
    id2seq = dict()

    seqfiles = []
    for ro, dirs, files in os.walk(seqdir, True):
        for i in files:
            if i.endswith('_segment.fa'):
                seqfiles.append(seqdir + i)
    for seqfile in seqfiles:
        f = open(seqfile, 'r')
        line = f.readline()
        id = ''
        while line:
            if line.startswith('>'):
                id = line.strip()
                id2seq[id] = ''
                line = f.readline()
                continue
            id2seq[id] += line.strip()
            line = f.readline()
        f.close()

    return id2seq

def extract_dna_outer_dna_inner(pos,searchlist,id2seq,framesize):
    """
    Extracting junction sequences by position
    :param pos: AS position list e.g. [3401, 4025, 4160, 4408, 4516, 5122, 5221, 5288, 5384, 5446, 5531, 5689]
    :param searchlist: targeted genomic segment sequences decided by boundary of AS
    :param id2seq: genomic segmented sequences
    :param framesize: chosen AS framesize, default 30bp*2
    :return: dna_out : list of outerside sequences
             dna_in : list of inside sequences
    """
    sequences = ''
    # print(searchlist)
    base_position = int(searchlist[0].strip().split('\t')[-1])

    for i in searchlist:
        # print(i)
        sequences += id2seq[i]

    dna_out = []
    dna_in = []
    # print(pos)
    # print(len(sequences))
    for i in range(len(pos)): # pos[1] number checked, minimal 51
        if i == 0: # first position
            if pos[i+1] - pos[i] + 1 < framesize : # RNA_interval_less_than_framesize

                if pos[i] < framesize: # DNA_head_position_less_than_framesize
                    dna_out.append('N' * (framesize - len(sequences[:pos[i]])) + sequences[:pos[i]]) # dna_left_outer
                    dna_in.append(sequences[pos[i]:pos[i+1]+1] + 'N'*(framesize-len(sequences[pos[i]:pos[i+1]+1])))  # dna_left_inner
                else: # DNA_head_position_larger_than_framesize
                    dna_out.append(sequences[pos[i]-base_position-framesize:pos[i]-base_position])
                    dna_in.append(sequences[pos[i]-base_position:pos[i + 1] - base_position + 1] + 'N'*(framesize-len(sequences[pos[i]-base_position:pos[i + 1] - base_position + 1])))

            else: # RNA_interval_larger_than_framesize

                if pos[i] < framesize:
                    dna_out.append('N' * (framesize - len(sequences[:pos[i]])) + sequences[:pos[i]])
                    dna_in.append(sequences[pos[i]:pos[i] + framesize])
                else:
                    dna_out.append(sequences[pos[i]-base_position-framesize:pos[i]-base_position])
                    dna_in.append(sequences[pos[i]-base_position:pos[i]-base_position + framesize])

        if i % 2 == 0 and i != 0: # even position
            if pos[i+1] - pos[i] + 1 < framesize : # RNA_interval_less_than_framesize

                if pos[i] - pos[i-1] < framesize: # DNA_interval_less_than_framesize
                    dna_out.append('N' * (framesize - len(sequences[pos[i-1]-base_position:pos[i]-base_position])) + sequences[pos[i-1]-base_position:pos[i]-base_position]) # dna_left_outer
                    dna_in.append(sequences[pos[i]-base_position:pos[i+1]-base_position+1] + 'N'*(framesize-len(sequences[pos[i]-base_position:pos[i+1]-base_position+1])))  # dna_left_inner
                else: # DNA_interval_larger_than_framesize
                    dna_out.append(sequences[pos[i]-base_position-framesize:pos[i]-base_position])
                    dna_in.append(sequences[pos[i]-base_position:pos[i + 1] - base_position + 1] + 'N'*(framesize-len(sequences[pos[i]-base_position:pos[i + 1] - base_position + 1])))

            else: # RNA_interval_larger_than_framesize

                if pos[i] - pos[i - 1] < framesize: # DNA_interval_less_than_framesize
                    dna_out.append('N' * (framesize - len(sequences[pos[i-1]-base_position:pos[i]-base_position])) + sequences[pos[i-1]-base_position:pos[i]-base_position]) # dna_left_outer
                    dna_in.append(sequences[pos[i]-base_position:pos[i]-base_position + framesize])
                else:# DNA_interval_larger_than_framesize
                    dna_out.append(sequences[pos[i]-base_position-framesize:pos[i]-base_position])
                    dna_in.append(sequences[pos[i]-base_position:pos[i]-base_position + framesize])

        if i % 2 == 1 and i != len(pos)-1: # odd position
            if pos[i]-pos[i-1]<framesize:# RNA_interval_less_than_framesize

                if pos[i+1]-pos[i]<framesize:# DNA_interval_less_than_framesize
                    dna_in.append('N'*(framesize-len(sequences[pos[i-1]-base_position:pos[i]-base_position+1]))+sequences[pos[i-1]-base_position:pos[i]-base_position+1])
                    dna_out.append(sequences[pos[i]-base_position+1:pos[i+1]-base_position]+'N'*(framesize-len(sequences[pos[i]-base_position+1:pos[i+1]-base_position])))
                else:# DNA_interval_larger_than_framesize
                    dna_in.append('N'*(framesize-len(sequences[pos[i-1]-base_position:pos[i]-base_position+1]))+sequences[pos[i-1]-base_position:pos[i]-base_position+1])
                    dna_out.append(sequences[pos[i]-base_position+1:pos[i]-base_position+framesize+1])

            else:# RNA_interval_larger_than_framesize

                if pos[i+1]-pos[i]<framesize:# DNA_interval_less_than_framesize
                    dna_in.append(sequences[pos[i]-base_position+1-framesize:pos[i]-base_position+1])
                    dna_out.append(sequences[pos[i]-base_position+1:pos[i+1]-base_position]+'N'*(framesize-len(sequences[pos[i]-base_position+1:pos[i+1]-base_position])))
                else:# DNA_interval_larger_than_framesize
                    dna_in.append(sequences[pos[i]-base_position+1-framesize:pos[i]-base_position+1])
                    dna_out.append(sequences[pos[i]-base_position+1:pos[i]-base_position+framesize+1])

        if i == len(pos)-1: # last position
            if pos[i]-pos[i-1]<framesize:# RNA_interval_less_than_framesize
                dna_in.append('N'*(framesize-len(sequences[pos[i-1]-base_position:pos[i]-base_position+1]))+sequences[pos[i-1]-base_position:pos[i]-base_position+1])
                dna_out.append(sequences[pos[i]-base_position+1:pos[i]-base_position+framesize+1])
            else:# RNA_interval_larger_than_framesize
                dna_in.append(sequences[pos[i]-base_position-framesize+1:pos[i]-base_position+1])
                dna_out.append(sequences[pos[i]-base_position+1:pos[i]-base_position+framesize+1])

    return dna_out,dna_in

def search_coordinator(position,framesize):
    # building sequences searching index
    search_dict = dict()
    for item in position:
        genome_id = '_'.join(item.strip().replace('|','_').split('_')[1:3])
        pos = position[item][-1]
        if pos[0] - framesize > 0:
            pos_min = (pos[0] - framesize)//8000
        else:
            pos_min = 0
        pos_max = (pos[-1] + framesize)//8000
        # print(item)
        # print(genome_id)
        # print(pos)
        # print(pos_min)
        # print(pos_max)
        # print('\n')
        search_id = []
        if pos_max == pos_min:
            search_id.append('>' + genome_id + '\t' + str(pos_min * 8000 + 1))  # >NC_001320.1	1
        else:
            for i in range(pos_max-pos_min + 1):
                search_id.append('>' + genome_id + '\t' + str((pos_min + i)*8000 + 1)) # >NC_001320.1	1
        search_dict[item] = search_id

    fw = open('./search_index_file.txt','w')
    for i in search_dict:
        fw.write(i+'\n')
        fw.write('\t'.join(search_dict[i])+'\n')
    fw.close()
    return search_dict

def load_junction_seqs(seqdir,positionfile,framesize):
    # output junction sequences
    position = load_position_info(positionfile)

    search_index = search_coordinator(position, framesize)

    id2seq = load_sequences_segment(seqdir)

    fw = open('./search_index_file_sequences.txt', 'w')
    for item in search_index:
        fw.write(item + '\n')
        id = ''
        seqs = ''
        for i in search_index[item]:
            id += i + '||'
            seqs += id2seq[i]
        fw.write(id + '\n')
        fw.write(seqs + '\n')
    fw.close()

    fw = open('./GCF_001433935.1_IRGSP-1.0_genomic.fna_DNA_outer_DNA_inner.fa', 'w')
    for item in position:
        if item in search_index:
            searchlist = search_index[item]
            # print(searchlist)
        else:
            continue
        comp = position[item][0]
        pos = position[item][1]
        dna_outer,dna_inner = extract_dna_outer_dna_inner(pos, searchlist, id2seq, framesize)
        # print(item)
        # print(position[item])
        # print(dna_outer)
        fw.write(item + '\t'+ comp + '\n')
        if comp == '1':
            re_inner = [reverse_seq(x) for x in dna_inner[::-1]]
            re_outer = [reverse_seq(x) for x in dna_outer[::-1]]
            for i in range(len(re_outer) // 2):
                fw.write(re_outer[2 * i]+'\t'+re_inner[2 * i]+'\t'+re_inner[2 * i + 1]+'\t'+re_outer[2 * i + 1]+'\n')
        else:
            for i in range(len(dna_outer) // 2):
                fw.write(dna_outer[2 * i]+'\t'+dna_inner[2 * i]+'\t'+dna_inner[2 * i + 1]+'\t'+dna_outer[2 * i + 1]+'\n')
    fw.close()


# def load_dna_outer(seqdir,positionfile,framesize):
#     position = dict()
#     id2seq = dict()
#
#     f = open(positionfile, 'r')
#     line = f.readline()
#     while line:
#         if line.startswith('>'):
#             id = line.strip().split('\t')[0]
#             comp = line.strip().split('\t')[1]
#             tmp = [int(x) for x in f.readline().strip().replace(',', '..').split('..')]
#             position[id] = [[comp],tmp]
#             # print(id)
#             # print(position[id])
#             # #  > lcl | NW_015379189.1_mrna_XM_015765349.1_50431
#             # [['0'], [3401, 4025, 4160, 4408, 4516, 5122, 5221, 5288, 5384, 5446, 5531, 5689]]
#
#         line = f.readline()
#     f.close()
#
#     search_index = search_coordinator(position,framesize)
#
#     seqfiles = []
#     for ro, dirs, files in os.walk(seqdir, True):
#         for i in files:
#             if i.endswith('_segment.fa'):
#                 seqfiles.append(seqdir + i)
#     for seqfile in seqfiles:
#         f = open(seqfile,'r')
#         line = f.readline()
#         id = ''
#         while line:
#             if line.startswith('>'):
#                 id = line.strip()
#                 id2seq[id] = ''
#                 line = f.readline()
#                 continue
#             id2seq[id] += line.strip()
#             line = f.readline()
#             # >NC_029260.1   1
#             # cctaaaccctaaaccctaaaccctaaaccctaaacc
#         f.close()
#
#     fw = open('./search_index_file_sequences.txt','w')
#     for item in search_index:
#         fw.write(item + '\n')
#         id = ''
#         seqs = ''
#         for i in search_index[item]:
#             id += i + '||'
#             seqs += id2seq[i]
#         fw.write(id + '\n')
#         fw.write(seqs + '\n')
#     fw.close()
#
#     fw = open('./GCF_001433935.1_IRGSP-1.0_genomic.fna_DNA_outer.fa', 'w')
#     for item in position:
#         if item in search_index:
#             searchlist = search_index[item]
#             print(searchlist)
#         else:
#             continue
#         # print(searchlist)
#         comp = position[item][0]
#         pos = position[item][1]
#         dna_outer = extract_dna_outer(pos,searchlist,id2seq,framesize)
#         # print(item)
#         # print(position[item])
#         # print(dna_outer)
#         fw.write(item + '\n')
#         if comp == '1':
#             temp = [x[::-1] for x in dna_outer[::-1]]
#             for i in range(len(temp)//2):
#                 fw.write(temp[2*i] + '\t' + temp[2*i+1] + '\n' )
#         else:
#             temp = dna_outer
#             for i in range(len(temp) // 2):
#                 fw.write(temp[2 * i] + '\t' + temp[2 * i + 1] + '\n')
#     fw.close()
#
#
# def extract_dna_outer(pos,searchlist,id2seq,framesize):
#     base_position = 0
#     sequences = ''
#     # print(searchlist)
#     base_position = int(searchlist[0].strip().split('\t')[-1])
#
#     for i in searchlist:
#         # print(i)
#         sequences += id2seq[i]
#
#     dna_out = []
#     # print(pos)
#     # print(len(sequences))
#     for i in range(len(pos)):
#         if i % 2 == 0 :
#             # left outer
#             if pos[i] < framesize:
#                 dna_out.append('N'*(framesize-len(sequences[:pos[i]])) + sequences[:pos[i]])
#             else:
#                 dna_out.append(sequences[pos[i]-base_position-framesize:pos[i]-base_position])
#         else:
#             # right outer
#             dna_out.append(sequences[pos[i]-base_position+1:pos[i]-base_position+framesize+2])
#
#     return dna_out






# def genomic_interval(file):
#     f = open(file,'r')
#     lines = f.readlines()
#     f.close()
#     interval = 100000
#     cnt = 0
#     for line in lines:
#         if line.startswith('>'):
#             pass
#         else:
#             tmp = [int(x) for x in line.strip().replace(',','..').split('..')]
#             for i in range(int(len(tmp)/2)):
#                 if i == 0:
#                     continue
#                 inter = tmp[2*i]-tmp[2*i-1] + 1
#                 if inter < 30:
#                     cnt += 1
#                 if interval > inter:
#                     interval = inter
#     print(interval)
#     print(cnt)
#     #  minimal interval = 3
#     #  for GCF_001433935.1_IRGSP-1.0_rna_from_genomic.fa_position.txt


def cut_genomic_file(file):
    # cut whole genome file into single fasta marker file
    f = open(file,'r')
    lines = f.readlines()
    f.close()

    seqs = {}
    id = ''
    for line in lines:
        if line.startswith('>'):
            id = line.strip().split(' ')[0]
            seqs[id] = []
            continue
        seqs[id].append(line)

    for item in seqs:
        fw = open('./genome_db/'+item.replace('>','')+'.fa','w')
        fw.write(item+'\n')
        for i in seqs[item]:
            fw.write(i)
        fw.close()

def load_db_list(root):
    # load single fasta files path
    filelist = []
    for ro, dirs, files in os.walk(root, True):
        for i in files:
            if i.startswith('N'):
                filelist.append(root + i)
    return filelist

def genome_coordinator(filelist):
    # cut file into every 8000bp
    for file in filelist:
        f = open(file, 'r')
        line = f.readline()

        fw = open(file+'_genome_segment.fa', 'w')
        id = ''
        cntline = 0
        lb = 0
        rb = 0
        seq = ''
        while line:
            if line.startswith('>'):
                cntline = 0
                lb = 0
                rb = 0
                id = line.strip().split(' ')[0]  # >NC_029256.1
            if cntline % 100 == 0 and cntline != 0:
                # fw.write(id + '\t' + str(lb + 1) + '_' + str(rb) + '\n')
                fw.write(id + '\t' + str(lb + 1) + '\n')
                fw.write(seq + '\n')
                lb = rb
                seq = ''
            line = f.readline()
            cntline += 1
            seq += line.strip()
            rb += len(line.strip())

        # fw.write(id + '\t' + str(lb + 1) + '_' + str(rb) + '\n')
        fw.write(id + '\t' + str(lb + 1) + '\n')
        fw.write(seq + '\n')

        fw.close()
        f.close()


def scan(file):
    f = open(file,'r')
    line = f.readline()
    max = 0
    min = 10000000
    while line:
        if line.startswith('>'):
            line = ' '.join([' '.join(x.split('..')) for x in f.readline().split('\t')[0].split(',')])
            num1 = int(line.split(' ')[0])
            num2 = int(line.split(' ')[-1])
            if num2 > num1:
                if num2 > max:
                    max = num2
                if num1 < min:
                    min = num1
            else:
                print('oo')
                if num1 > max:
                    max = num1
                if num2 < min:
                    min = num2
        line = f.readline()
    print(max)
    print(min)

def count_position_number(file):
    f = open(file,'r')
    lines = f.readlines()
    f.close()

    for line in lines:
        if line.startswith('>'):
            continue
        else:
            tmp = line.strip().replace(',','..').split('..')
            if len(tmp) % 2 ==0:
                if int(tmp[1])<100:
                    print(line)
            else:
                print(line)

def negative_sample(seqdir,positionfile,framesize):
    position = load_position_info(positionfile)

    search_index = search_coordinator(position, framesize)

    id2seq = load_sequences_segment(seqdir)

    AS_segment = set()
    genomic_segment = set()
    for item in search_index:
        for j in search_index[item]:
            AS_segment.add(j)
    for item in id2seq:
        genomic_segment.add(item)

    # print(len(genomic_segment)) # 46833
    # print(len(AS_segment)) # 28338
    # print(len(genomic_segment-AS_segment)) # 18495
    temp = dict()
    for seqid in genomic_segment-AS_segment:
        nega_seq = cut_negative_seq(id2seq[seqid],framesize)
        temp[seqid] = nega_seq

    fw= open('./negative_sample.txt','w')
    for item in temp:
        for i in range(len(temp[item])//2):
            fw.write(item+'\t'+'junction_'+str(i+1)+'\t')
            fw.write(temp[item][2*i]+'\t'+temp[item][2*i+1]+'\t'+'0'+'\n')
    fw.close()


def cut_negative_seq(seq,framesize):
    index = np.random.choice(len(seq)-2*framesize-1,56,replace=False)
    nega_seq = []
    for i in index:
        nega_seq.append(seq[i:i+framesize]+'\t'+seq[i+framesize:i+2*framesize])
    return nega_seq





def chech_title(file):
    f = open(file, 'r')
    lines = f.readlines()
    f.close()
    header = set()
    for line in lines:
        if line.startswith('>'):
            header.add(line.strip().split('_')[2])
    print(header)
    # {'mrna', 'miscrna', 'ncrna', 'trna'}

def positive_sample(file):
    f = open(file, 'r')
    line = f.readline()

    AS_pool = dict()
    id = ''
    while line:
        if line.startswith('>'):
            id = line.strip()
            AS_pool[id] = []
            line = f.readline()
        AS_pool[id].append(line.strip())
        line = f.readline()
    f.close()

    fw = open(file+'_positive_sample.txt', 'w')
    for item in AS_pool:
        if item.split('_')[2] == 'mrna':
            for i in range(len(AS_pool[item])):
                fw.write(item + '\t' + 'junction_'+str(i+1)+'\t')
                fw.write(AS_pool[item][i]+'\t'+'1'+'\n')
    fw.close()
    # output more than 257196 AS positive item


if __name__=='__main__':
    # load_position('/home/david/Desktop/asdecoder/GCF_001433935.1_IRGSP-1.0_rna_from_genomic.fa')
    # load_seq('/home/david/Desktop/asdecoder/GCF_001433935.1_IRGSP-1.0_rna_from_genomic.fa')

    # scan('/home/david/Desktop/asdecoder/GCF_001433935.1_IRGSP-1.0_rna_from_genomic.fa_position.txt')

    # load_rna_inner('/home/david/Desktop/asdecoder/GCF_001433935.1_IRGSP-1.0_rna_from_genomic.fa_seq.txt',
    #                '/home/david/Desktop/asdecoder/GCF_001433935.1_IRGSP-1.0_rna_from_genomic.fa_position.txt',
    #                30)

    # genomic_interval('./GCF_001433935.1_IRGSP-1.0_rna_from_genomic.fa_position.txt')
    # genome_coordinator('/home/david/Desktop/asdecoder/genome_db/GCF_001433935.1_IRGSP-1.0_genomic.fna')

    # load_dna_outer('/home/david/Desktop/asdecoder/genome_db/GCF_001433935.1_IRGSP-1.0_genomic.fna_genome_segment.fa',
    #                '/home/david/Desktop/asdecoder/GCF_001433935.1_IRGSP-1.0_rna_from_genomic.fa_position.txt',
    #                30)



    # load_junction_seqs('./genome_db/segment/',
    #                    './GCF_001433935.1_IRGSP-1.0_rna_from_genomic.fa_position.txt',
    #                    30)

    # count_position_number('./GCF_001433935.1_IRGSP-1.0_rna_from_genomic.fa_position.txt')

    # Cut genomic file and build coordinator
    # cut_genomic_file('./genome_db/GCF_001433935.1_IRGSP-1.0_genomic.fna')
    # filelist = load_db_list('./genome_db/')
    # genome_coordinator(filelist)

    negative_sample('./genome_db/segment/',
                    './GCF_001433935.1_IRGSP-1.0_rna_from_genomic.fa_position.txt',
                    30)

    # chech_title('./GCF_001433935.1_IRGSP-1.0_genomic.fna_DNA_outer_DNA_inner.fa')

    # positive_sample('./GCF_001433935.1_IRGSP-1.0_genomic.fna_DNA_outer_DNA_inner.fa')