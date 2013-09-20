#-- Utilities
#
# Some general-purpose tools
# Author: Jason Winget
# Version: 0.1
#
#--

import time

def Timer(fn):
	''' Decorator to time a function's execution '''
	def timed(*args, **kw):
		start = time.time()
		result = fn(*args, **kw)
		elapsed = str(round(time.time() - start, 2))
		print fn.__name__+' completed in '+elapsed+' s'
		return result
	return timed

# Binary seach tree implementation
# From http://www.laurentluce.com/posts/binary-search-tree-library-in-python/

class Node:
	''' Tree node: left and right child + data '''
	def __init__(self, data):
		self.left = None
		self.right = None
		self.data = data
	def insert(self, data):
		''' Recursive function to insert a new node in the tree '''
		if data < self.data:
			if self.left is None:
				self.left = Node(data)
			else:
				self.left.insert(data)
		else:
			if self.right is None:
				self.right = Node(data)
			else:
				self.right.insert(data)
	def lookup(self, data, parent = None):
		''' Lookup a node in the tree '''
		if data < self.data:
			if self.left is None:
				return None, None
			return self.left.lookup(data, self)
		elif data > self.data:
			if self.right is None:
				return None, None
			return self.right.lookup(data, self)
		else:
			return self, parent
	def children_count(self):
		''' Return the number of children of a node '''
		if node is None:
			return None
		cnt = 0
		if self.left:
			cnt += 1
		if self.right:
			cnt += 1
		return cnt
	def delete(self, data):
		''' Remove a node from the tree '''
		node, parent = self.lookup(data)
		if node is not None:
			children_count = node.children_count()
			if children_count == 0:
				if parent.left is node:
					parent.left = None
				else:
					parent.right = None
				del node
			elif chidren_count == 1:
				if node.left:
					n = node.left
				else:
					n = node.right
				if parent:
					if parent.left is node:
						parent.left = n
					else:
						parent.right = n
				del node
			else:
				parent = node
				successor = node.right
				while successor.left:
					parent = successor
					successor = successor.left
				node.data = successor.data
				if parent.left == successor:
					parent.left = successor.right
				else:
					parent.right = successor.right
	def print_tree(self):
		''' Print the tree content '''
		if self.left:
			self.left.print_tree()
		print self.data,
		if self.right:
			self.right.print_tree()
