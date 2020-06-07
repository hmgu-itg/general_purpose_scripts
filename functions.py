import re

# this function is primarily used to get variant's position
def parseSPDI(string):
    L=string.rsplit(":")
    c=L[0]
    m=re.search("NC_0+(\d+)\.\d+",L[0])
    if m:
        c=m.group(1)
    pos=int(L[1])
    ref=L[2]
    alt=L[3]
    lref=len(ref)
    lalt=len(alt)

    # ref can be the length of the deleted sequence
    m=re.search("^(\d+)$",ref)
    if m:
        if lalt==1: # SNP
            pos=pos+1
    elif lref==1 and lalt==1:# SNP
        pos=pos+1

    return {"chr":c,"pos":pos,"ref":ref,"alt":alt}
