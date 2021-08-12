from collections import  Counter

class Node:
    def __init__(self, char=""):
        self.right = None
        self.left = None
        self.char = char

class Queue:
    def __init__(self):
        self.queue = []
        self.size = 0

    def get_item(self):
        self.size -= 1
        return self.queue.pop(0)

    def input_item(self, item, priority):
        if not self.queue:
            self.queue = [(item, priority)]
            self.size += 1
        else:
            for i in range(len(self.queue)):
                if priority <= self.queue[i][1]:
                    self.queue.insert(i, (item,priority))
                    self.size += 1
                    return
            self.queue.append((item, priority))
            self.size += 1

def buildTree(string):
    queue = Queue()
    # 出现的字符按照ASCII码进行排序
    probability = {char: priority for (char, priority) in (sorted(Counter(string).items(), key=lambda item:item[0]))}

    for char, priority in probability.items():
        node = Node(char)
        queue.input_item(node, priority)

    while queue.size:
        if queue.size == 1:
            return queue.get_item()[0]
        temp = Node()
        node1 = queue.get_item()
        node2 = queue.get_item()
        temp.left = node1[0]
        temp.right = node2[0]
        queue.input_item(temp, node1[1]+node2[1])

def traversal(node, table, code=""):
    if not node.left and not node.right:
        table[node.char] = code
    if node.left:

        traversal(node.left,table, code+"0")
    if node.right:
        traversal(node.right,table, code+"1")

def buildTable(root):
    table = {}
    traversal(root, table)
    return table


def encode(string, table):
    encoded = ""
    for c in string:
        encoded += table[c]

    return bytes(int(x) for x in encoded)


def decode(string, tree):
    curse = tree
    decoded = ""
    for c in string:
        if not curse.left and not curse.right:
            decoded += curse.char
            curse = tree
        if c == '1':
            curse = curse.right
        if c == '0':
            curse = curse.left
    return decoded

def str2byte(str):
    return bytes(int(c) for c in str)

def byte2str(byte):
    return ''.join(str(b)for b in byte)


if __name__ == '__main__':
    import sys
    s1 = "I love FishC.com!"
    s2 = "0011111000111"
    root = buildTree(s1)
    table = buildTable(root)
    encoded = encode(s1, table)
    decoded = decode(s2, root)
    print(encoded == str2byte('100011001110010000101011010010100000101011111011111101100101101110'))
    print(decoded == "oCo")
    print(sys.getsizeof(str2byte('100011001110010000101011010010100000101011111011111101100101101110')))




