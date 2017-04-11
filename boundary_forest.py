import skbio
from numpy import *
from numpy.random import *
import Node
from SW_distance import *
import datetime

begin=datetime.datetime.now()

#read fasta into a dictionary
def read_sequences(filename):
    f=open(filename)
    lines=f.readlines()
    f.close()
    titles={}
    mytitle=""
    for i in lines:
        if i[0]==">":
            mytitle=i.strip()
            titles[mytitle]=""
        else:
            titles[mytitle]+=(i.upper()).strip()
    return titles

def NW_distance(prot1ID,prot2ID):
    prot1=skbio.sequence.Protein(prot1ID)
    prot2=skbio.sequence.Protein(prot2ID)
    return skbio.alignment.global_pairwise_align_protein(prot1,prot2)[1]

def dict_to_list(mydictionary):
    mylist=[[] for i in mydictionary]
    j=0
    for i in mydictionary:
        mylist[j]=[i,mydictionary[i]]
        j+=1
    return mylist

def norm(v1,v2):
    if len(v1)!=len(v2):
        print("vectors need to have the same length")
        return
    vdif=v1
    for i in range(0,len(v2)):
        vdif[i]-=v2[i]
    return (sum([i**2 for i in vdif]))**0.5

def SW_distance_weighed(seq1,seq2):
    sw_raw=SW_distance(seq1,seq2)
    weight=min(len(seq1),len(seq2))
    return sw_raw/weight

def length_match(seq1,seq2):
    l1=len(seq1)
    l2=len(seq2)
    if 1.0*min(l1,l2)/max(l1,l2)<0.7: return False
    else: return True
#get possible children based on length of sequences
def get_possible_children(seq1, parentnode):
    children=parentnode.children
    possible_children=[]
    for child in children:
        childseq=data[child.data][1]
        if length_match(seq1,childseq):
            possible_children.append(child)
    if possible_children==[]: return children
    return possible_children

def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""
#pull out the protein annotation. Return true only if the two are an exact match
#and NOT hypothetical.protein
def prot_match(title1,title2):
    prot1=find_between(title1,"[protein=","]")
    prot2=find_between(title2,"[protein=","]")
    if prot1=="hypothetical.protein":
        return False
    else:
        return prot1==prot2
    
#threshold for comparison of label vectors
eps=0
#distance metric for comparing positions
dis=NW_distance
#max # children in the tree
max_deg=10

n=100
p=10

#get real sequence data
seqs="C:/Users/defne/Documents/BC/TVO/pangenome/clustering_test_proteins.fasta"
data=read_sequences(seqs)
data=dict_to_list(data)
n=len(data)

#data=rand(n,p)
data_order_ix=range(0,n)

tree=Node.Node(data_order_ix[0],0)
secondnode=Node.Node(data_order_ix[1],1)
tree.add_child(secondnode)
num_tree_nodes=2

for i in range(2,n):
    ele_being_proc_data_ix=data_order_ix[i]
    current_tree_node=tree
    while True:
        curr_node_data_ix=current_tree_node.data
        curr_node_all_children_node_ix=current_tree_node.children
        print(ele_being_proc_data_ix,data[ele_being_proc_data_ix][0],data[curr_node_data_ix][0])
        ##find the child closes to query
        smallest_dist=10000000
        best_child_node=None
        
        #compare to parent by name and by length
        v1=data[ele_being_proc_data_ix][1]
        v1t=data[ele_being_proc_data_ix][0]
        v3=data[curr_node_data_ix][1]
        v3t=data[curr_node_data_ix][0]
        lengths=length_match(v1,v3)
        same_protein=prot_match(v1t,v3t)
        
        #don't even bother if they have the same annotation
        if same_protein:break
        #if there are no children, add the element as a child to this node
        elif len(current_tree_node.children)==0:
            num_tree_nodes+=1
            ele_being_proc_node=Node.Node(ele_being_proc_data_ix,num_tree_nodes)
            current_tree_node.add_child(ele_being_proc_node)
            break

        #find the best child
        for child in  curr_node_all_children_node_ix:
            ##distance function - compare to children
            v1=data[ele_being_proc_data_ix][1]
            v2=data[child.data][1]
            dis=SW_distance_weighed(v1,v2)
            if dis<smallest_dist:
                smallest_dist=dis
                best_child_node=child
        #/for
                
        #check if lengths are similar and whether curent node has space for children
        if lengths and len(current_tree_node.children)<max_deg:
            ##only compute SW if current node has space for a new
            ##child AND if the sequences are similar in length
            dis2=SW_distance_weighed(v1,v3)
            ##ignore the point if it's within a threshold distance
            if dis2<smallest_dist and dis2>eps:
                num_tree_nodes+=1
                ele_being_proc_node=Node.Node(ele_being_proc_data_ix,num_tree_nodes)
                current_tree_node.add_child(ele_being_proc_node)
            #/end if
            break
        else:
            current_tree_node=best_child_node
    ##/while
##/for

ending=datetime.datetime.now()

print((ending-begin)/datetime.timedelta(minutes=1))

tree.print_tree()

##pairwise=[[0 for i in range(0,25)] for i in range(0,25)]
##for i in range(0,25):
##    for j in range(0,25):
##        pairwise[i][j]=SW_distance_weighed(data[i][1],data[j][1])
##
##for i in pairwise: print (i)
