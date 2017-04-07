import skbio
from numpy import *
from numpy.random import *
import Node
from SW_distance import *

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
	
#threshold for comparison of label vectors
eps=0.1
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
        print(ele_being_proc_data_ix,data[ele_being_proc_data_ix][0])
        curr_node_data_ix=current_tree_node.data
        curr_node_all_children_node_ix=current_tree_node.children
            
        ##find the child closes to query
        smallest_dist=10000000
        best_child_node=None
        for child in  curr_node_all_children_node_ix:
            ##distance function - compare to children
            v1=data[ele_being_proc_data_ix][1]
            v2=data[child.data][1]
            #dis=NW_distance(v1,v2)
            dis=SW_distance(v1,v2)
            if dis<smallest_dist:
                smallest_dist=dis
                best_child_node=child
        #/for
            
        #distance function - compare to parent
        v1=data[ele_being_proc_data_ix][1]
        v3=data[curr_node_data_ix][1]
        #dis2=NW_distance(v1,v3)
        dis2=SW_distance(v1,v3)
        
        #compare dis to dis2. whoever is smaller gets element 
        #being processed as a child
        if dis2<smallest_dist and len(current_tree_node.children)<max_deg:
            ##ignore the point if it's within a threshold distance
            if dis2>eps:
                num_tree_nodes+=1
                ele_being_proc_node=Node.Node(ele_being_proc_data_ix,num_tree_nodes)
                current_tree_node.add_child(ele_being_proc_node)
            #/end if
            break
        else:
            current_tree_node=best_child_node
    ##/while
##/for

tree.print_tree()
