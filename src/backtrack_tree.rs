use slab::Slab;

#[derive(Copy, Clone, Debug)]
pub struct NodeId(usize);

/// Very simply tree structure. Except for root nodes, each node has a parent.
/// The only direction of traversal is from children to parents (towards a root).
/// Other relationships are not stored. It's possible to build and traverse trees with multiple roots.
#[derive(Default)]
pub struct Tree<T>(Slab<Node<T>>);

impl<T> Tree<T> {
    /// Create a new tree
    pub fn new() -> Self {
        Self(Slab::new())
    }

    /// Create a new tree with specified capacity to avoid frequent allocations
    /// if the maximal size of the tree is known upfront
    pub fn with_capacity(capacity: usize) -> Self {
        Self(Slab::with_capacity(capacity))
    }

    /// Remove node from the tree. There is no safety net, it's totally possible to remove the root of a subtree without any warning.
    /// Also, the key of the removed item will be reused. So please ensure to invalidate the old key.
    pub fn remove(&mut self, key: NodeId) {
        let _ = self.0.remove(key.0);
    }

    /// Add node to tree under `parent`. If `parent` is `None`, the current node is inserted as a root.
    /// It's possible to insert multiple roots into one tree.
    pub fn add_node<U>(&mut self, value: T, parent: U) -> NodeId
    where
        U: Into<Option<NodeId>>,
    {
        NodeId(self.0.insert(Node {
            value,
            parent: parent.into(),
        }))
    }

    /// Returns an _inclusive_ iterator over the data associated with the ancestors of a given node ID
    pub fn ancestors(&self, id: NodeId) -> Ancestors<T> {
        Ancestors {
            tree: self,
            state: Some(id),
        }
    }

    /// Returns the number of elements in the tree
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns if the tree is empty
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the number of elements the Tree can hold without reallocating
    pub fn capacity(&self) -> usize {
        self.0.capacity()
    }

    /// Clears the tree, keeps the allocated memory for reuse
    pub fn clear(&mut self) {
        self.0.clear()
    }
}

/// The node struct is not exposed at the API surface
struct Node<T> {
    value: T,
    parent: Option<NodeId>,
}

/// Inclusive iterator over the data associated with the ancestors of a given node ID
pub struct Ancestors<'a, T> {
    tree: &'a Tree<T>,
    state: Option<NodeId>,
}

impl<'a, T> Iterator for Ancestors<'a, T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(state) = &self.state {
            if let Some(node) = self.tree.0.get(state.0) {
                self.state = node.parent;
                return Some(&node.value);
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tree() {
        let value = 15;
        let mut tree = Tree::new();

        let root_id = tree.add_node(value, None);
        let second_id = tree.add_node(value + 1, root_id);
        let third_id = tree.add_node(value + 2, second_id);
        let fourth_id = tree.add_node(value + 3, third_id);

        let out = tree.ancestors(fourth_id).copied().collect::<Vec<_>>();
        assert_eq!(out, vec![18, 17, 16, 15]);

        tree.remove(second_id);
        let out = tree.ancestors(fourth_id).copied().collect::<Vec<_>>();
        assert_eq!(out, vec![18, 17]);
    }

    #[test]
    fn test_length() {
        let value = 15;
        let mut tree = Tree::new();
        tree.add_node(value, None);
        assert_eq!(tree.len(), 1);
        assert_eq!(tree.len(), tree.0.len());
    }

    #[test]
    fn test_capacity() {
        let value = 15;
        let mut tree = Tree::with_capacity(5);
        tree.add_node(value, None);
        assert_eq!(tree.capacity(), 5);
        assert_eq!(tree.capacity(), tree.0.capacity());
    }

    #[test]
    fn test_clear() {
        let value = 15;
        let mut tree = Tree::new();
        let root_id = tree.add_node(value, None);
        tree.add_node(value, root_id);
        assert_eq!(tree.len(), 2);
        tree.clear();
        assert_eq!(tree.len(), 0);
    }
}
