class Node(object):
    "Generic tree node."
    def __init__(self, data, index, level=0, parent=None, children=None):
        self.data = data
        self.index = index
        self.parent = parent
        self.level = level
        self.children = []
        if children is not None:
            for child in children:
                self.add_child(child)
    def add_child(self, node):
        self.children.append(node)
        node.parent=self
        node.level=self.level+1
    def is_leaf(self):
        return self.children==[]
    def print_node(self,data=False):
        indent=self.level*" "
        if data==False:
            nodestr=str(self.index)+"\n"
        else:
            nodestr=str(self.data)+"\n"
        print(indent+nodestr)
    def print_tree(self):
        current_level=self.level
        current_node=self
        if current_node.is_leaf()==False:
            current_node.print_node()
            for childnode in current_node.children:
                childnode.print_tree()
        else:
            current_node.print_node()


x=Node(12,1)
y=Node(30,2,x)
x.add_child(y)
x.print_tree()
