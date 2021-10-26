'''
class for dbcan-profiler 

'''

### paf record sample
'''
1	string	Query sequence name
2	int	Query sequence length
3	int	Query start (0-based; BED-like; closed)
4	int	Query end (0-based; BED-like; open)
5	char	Relative strand: "+" or "-"
6	string	Target sequence name
7	int	Target sequence length
8	int	Target start on original strand (0-based)
9	int	Target end on original strand (0-based)
10	int	Number of residue matches
11	int	Alignment block length
12	int	Mapping quality (0-255; 255 for missing)
13      attribute

S0R223/1	149	0	149	+	AGR60979.1|GT8	1014	627	776	145	149	1	AS:i:-8	XS:i:-8	XN:i:0	XM:i:4	XO:i:0	XG:i:0	NM:i:4	MD:Z:121A11A7A5G1	YT:Z:UU	cg:Z:149M
S0R274/1	150	0	150	-	APC33536.1|GH42	2004	918	1068	148	150	1	AS:i:-4	XS:i:-4	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:11G20C117	YT:Z:UU	cg:Z:150M
S0R318/1	150	0	150	-	AGA56765.1|GH9	1863	667	817	147	150	42	AS:i:-6	XN:i:0	XM:i:3	XO:i:0	XG:i:0	NM:i:3	MD:Z:11G32G13C91	YT:Z:UU	cg:Z:150M
'''

### filter the CAZy 

def CAZy_filter(cazy):
    return set([aa.split("_")[0] for aa in cazy])

### 
class PafRecord(object):
    def __init__(self,lines):
        ### basic information 
        self.Qsn = lines[0]
        self.Qsl = lines[1]
        self.Qs  = lines[2]
        self.Qe  = lines[3]
        self.Strand = lines[4]
        self.Tsn = lines[5]
        self.Tsl = lines[6]
        self.Ts  = lines[7]
        self.Te  = lines[8]
        self.Nrm = lines[9]
        self.Abl = lines[10]
        self.Mq  = lines[11] ### if the paf was converted from sam, Mq here stands for the MAPQ 
        
        ### deal information
        self.SeqID = self.Tsn.split('|')[0]
        #self.CAZys = self.Tsn.strip("|").split("|")[1:] ### seqid|cazy1|cazy2|...| ## include subfamily
        self.CAZys = CAZy_filter(self.Tsn.strip("|").split("|")[1:]) ### seqid|cazy1|cazy2|...| ## not subfamily
        self.UniReadId = lines[0].split("/")[0]

    def __str__(self):
        return "\t".join([getattr(self, value) for value in vars(self) if value != "CAZys"])

### design for assemble free method to deal with the paf format file
### ### paf: https://github.com/lh3/miniasm/blob/master/PAF.md
### sam format file also convert to paf 
### https://bioconvert.readthedocs.io/en/master/developer_guide.html#how-to-update-bioconvert-on-bioconda
### diamond blastx: output format paf

class Paf(object):
    def __init__(self,filename):
        self.records = [PafRecord(line.split()) for line in open(filename)]
    def __iter__(self):
        return iter(self.records)
    ### get reads id
    def GetReadId(self):
        return [record.Qsn for record in self]
    ### get protein id
    def GetSeqId(self): 
        return [record.SeqID for record in self]
    ### get protein id: protein length dictory
    def GetSeqLen(self):
        return {record.SeqID:record.Tsl for record in self}
    ### get CAZy family id 2 protein id: one-many
    def CAZy2SeqID(self,CazySeqId):
        for record in self:
            for cazy in record.CAZys:
                CazySeqId.setdefault(cazy,[]).append(record.SeqID)
    ## get protein id 2 read is: one-many
    def SeqID2ReadID(self,aa):
        for record in self:
            aa.setdefault(record.SeqID,[]).append(record.Qsn)
    def ReadID2Record(self):
        return {record.Qsn:record for record in self}
    def Output(self):
        [print (record) for record in self]
    ## the CAZy information for megahit are not Qsn instead of they are in the 
    def Assign_CAZy_megahit(self):
        for cazy in self:
            cazy.CAZys = CAZy_filter(cazy.Qsn.strip("|").split("|")[1:])

class Golden_smaple_record(object):
    def __init__(self,lines):
        self.OTU = lines[0]
        self.genomeid = lines[1]
        self.seqid = lines[2]
    def __str__(self):
        return "\t".join([getattr(self, value) for value in vars(self)])

class Golden_smaple(object):
    def __init__(self,filename):
        self.records = [Golden_smaple_record(line.split()) for line in open(filename)]
    def __iter__(self):
        return iter(self.records)
    def OTUdict(self):
        dict1 = {}
        for record in self:
            dict1.setdefault(record.OTU,[]).append(record)
        return dict1

class Abundance_record(object):
    def __init__(self,lines):
        self.genomeid = lines[0]
        self.abund = lines[1]
    def __str__(self):
        return "\t".join([getattr(self, value) for value in vars(self)])

class Abundance(object):
    def __init__(self,filename):
        self.records = [Abundance_record(line.split()) for line in open(filename) if float(line.split()[1]) > 0]
    def __iter__(self):
        return iter(self.records)
    def __str__(self):
        return "\n".join([str(record) for record in self if record])
    def __len__(self):
        return len(self.records)


#taxid   kindom  phylum  class   order   family  genus   species
class Taxnomony_record(object):
    def __init__(self,lines):
        self.taxid,self.kindom,self.phylum,self.tclass,self.order,self.family,self.genus,self.species = lines
    def __str__(self):
        return ",".join([getattr(self, value) for value in vars(self)])

class Taxnomony(object):
    def __init__(self,filename):
        self.records = [Taxnomony_record(line.rstrip("\n").split("\t")) for line in open(filename)]
    def __iter__(self):
        return iter(self.records)
    def __str__(self):
        return "\n".join([str(record) for record in self if record])
    def __len__(self):
        return len(self.records)
    def taxid2tttt(self):
        return {record.taxid:record for record in self}
