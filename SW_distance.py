
#grab the blosum62 scores for amino acid comparison
blosum62file="C:/Users/defne/Documents/BC/TVO/pangenome/BF_Jose/blosum62.txt"
f=open(blosum62file)
lines=f.readlines()
f.close()
aas = lines[0].split()
blosumdict={}
#make blosum dictionary
for line in lines[1:]:
    blosumrow = line.split()
    rowdict={}
    AA=blosumrow[0]
    for i in range(1,len(blosumrow)):
        rowdict[aas[i-1]]=int(blosumrow[i])
    blosumdict[AA]=rowdict

#smith waterman distance - given seq1 and seq2 compute an alignment score
#gap opening and extending penalties are different
#similarity/match=lower score
def SW_distance(seq1, seq2,gap_open=10,gap_extend=0.5):
    rows=len(seq1)+1
    cols=len(seq2)+1

    #make the score and path matrices
    score_matrix=[[0 for col in range(0,cols)] for row in range(0,rows)]
    path_matrix=[["" for col in range(0,cols)] for row in range(0,rows)]
    #the first row and first column will simply be gap opening and extending
    for i in path_matrix[0][1:]:
        i="l"
    for row in path_matrix[1:]:
        row[0]="u"
    score_matrix[0][1]=gap_open
    score_matrix[1][0]=gap_open
    for i in range(2,cols):
        score_matrix[0][i]=score_matrix[0][i-1]+gap_extend
    for j in range(2,rows):
        score_matrix[j][0]=score_matrix[j-1][0]+gap_extend

    #keep track of maximum score, this is the metric to be returned
    max_score=0
    max_pos=(0,0)

    #update the scores one by one
    for i in range(1,rows):
        for j in range(1,cols):
            current_score,current_dir=calculate_score(score_matrix,path_matrix,i,j,seq1,seq2,gap_open,gap_extend)
            score_matrix[i][j]=current_score
            path_matrix[i][j]=current_dir
            if current_score>max_score:
                max_score=current_score
                max_pos=(i,j)
    return max_score

#calculate the score of the current position in the matrix
def calculate_score(mymatrix,mypath,row,col,seq1,seq2,gap_open,gap_extend):
    char1=seq1[row-1]
    char2=seq2[col-1]
    #make sure to use the appropriate gap penalty depending on 
    #whether we are opening or extending a gap
    if mypath[row-1][col]=="u":
        gap_penalty_up=gap_extend
    else:
        gap_penalty_up=gap_open
        
    if mypath[row][col-1]=="l":
        gap_penalty_left=gap_extend
    else:
        gap_penalty_left=gap_open

    #3 possible ways to get to the current position.
    #UP: gap on seq2
    #LEFT: gap on seq1
    #DIAGONAL: no gap, amino acid match
    score_from_up=mymatrix[row-1][col]+gap_penalty_up
    score_from_left=mymatrix[row][col-1]+gap_penalty_left
    score_from_diag=mymatrix[row-1][col-1]-match_score(char1,char2)
    
    thisscore= max(score_from_up,score_from_left,score_from_diag)
    thisdir=""
    
    #keep track of the direction of movement so we know whether
    #we open or extend gaps in the future.
    if thisscore==score_from_up:
        thisdir="u"
    elif thisscore==score_from_left:
        thisdir="l"
    else:
        thisdir="d"
    return thisscore, thisdir


def match_score(aa1, aa2):
    return blosumdict[aa1][aa2]
