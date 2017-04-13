from numpy import *
from numpy.random import *
import ArrayTree
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

def dict_to_list(mydictionary):
    mylist=[[] for i in mydictionary]
    j=0
    for i in mydictionary:
        mylist[j]=[i,mydictionary[i]]
        j+=1
    return mylist

def SW_distance_weighed(seq1,seq2):
    sw_raw=SW_distance(seq1,seq2)
    weight=min(len(seq1),len(seq2))
    return sw_raw/weight

def norm(v1,v2):
	if len(v1)!=len(v2):
		print("vectors need to have the same length")
		return
	vdif=v1
	for i in range(0,len(v2)):
		vdif[i]-=v2[i]
	return (sum([i**2 for i in vdif]))**0.5

def eucl(v1,v2):
    return((v1[0]-v2[0])**2 + (v1[1]-v2[1])**2)**0.5

def read_toy_data(file):
    f=open(file)
    lines=f.readlines()
    f.close()
    vectors=[]
    for i in lines:
        vectors.append([float(j) for j in i.split()])
    return vectors


#threshold for comparison of label vectors
eps=0
#max # children in the tree
max_deg=10

'''
data=read_toy_data("C:/Users/defne/Documents/BC/TVO/pangenome/BF_Jose/testdata.txt")
n=len(data)
data_order_ix=range(0,n)
'''

#get real sequence data
seqs="C:/Users/defne/Documents/BC/TVO/pangenome/clustering_test_proteins.fasta"
data=read_sequences(seqs)
data=dict_to_list(data)
n=len(data)
#data=rand(n,p)
data_order_ix=range(0,n)

tree=ArrayTree.ArrayTree(0,data_order_ix[0])
tree.add_child(0,[1,data_order_ix[1],[],0])
num_tree_nodes=2

for i in range(2,n):
        ele_being_proc_data_ix=data_order_ix[i]
        current_tree_node_ix=0

        while True:
                current_tree_node=tree.tree[current_tree_node_ix]
                curr_node_data_ix=current_tree_node[1]
                curr_node_all_children_node_ix=current_tree_node[2]
                print(ele_being_proc_data_ix,data[ele_being_proc_data_ix][0],data[curr_node_data_ix][0])
                ##find the child closes to query
                smallest_dist=10000000
                best_child_node=None
                for child in curr_node_all_children_node_ix:
                        ##distance function - compare to children
                        v1=data[ele_being_proc_data_ix][1]
                        v2=data[tree.tree[child][1]][1]
                        dis=SW_distance_weighed(v1,v2)
                        if dis<smallest_dist:
                                smallest_dist=dis
                                best_child_node=child
                #/for

                #distance function - compare to parent
                v1=data[ele_being_proc_data_ix][1]
                v3=data[current_tree_node[1]][1]
                dis2=SW_distance_weighed(v1,v3)

                #compare dis to dis2. whoever is smaller gets element 
                #being processed as a child
                if dis2<smallest_dist and len(current_tree_node[2])<max_deg:
                        ##ignore the point if it's within a threshold distance
                        if dis2>eps:
                                num_tree_nodes+=1
                                tree.add_child(current_tree_node_ix,[num_tree_nodes-1,ele_being_proc_data_ix,[],0])
                        #/end if
                        break
                else:
                        current_tree_node_ix=best_child_node
        ##/while
##/for

tree.print_tree()

ending=datetime.datetime.now()

print((ending-begin)/datetime.timedelta(minutes=1))
