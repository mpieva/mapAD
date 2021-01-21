use std::convert::TryFrom;

use slab::Slab;

/// Type wrapper around `u32` instead of `usize`
/// to save memory
#[derive(Copy, Clone, Debug)]
pub struct NodeId(u32);

/// Very simply tree structure. Except for root nodes, each node has a parent.
/// The only direction of traversal is from children to parents (towards a root).
/// Other relationships are not stored. It's possible to build and traverse trees with multiple roots.
#[derive(Default)]
pub struct Tree<T>(Slab<Node<T>>);

impl<T> Tree<T>
where
    T: Default,
{
    /// Create a new tree
    pub fn new() -> Self {
        Self(Slab::new())
    }

    /// `NodeId(0)` is reserved for the root node
    fn add_internal_root(&mut self) {
        let id = self.0.insert(Node {
            value: T::default(),
            parent: NodeId(0),
        });
        assert_eq!(id, 0);
    }

    /// Create a new tree with specified capacity to avoid frequent allocations
    /// if the maximal size of the tree is known upfront
    pub fn with_capacity(capacity: u32) -> Self {
        Self(Slab::with_capacity(capacity as usize))
    }

    /// Remove node from the tree. There is no safety net, it's totally possible to remove the root of a subtree without any warning.
    /// Also, the key of the removed item will be reused. So please ensure to invalidate the old key.
    pub fn remove(&mut self, key: NodeId) {
        if key.0 != 0 {
            let _ = self.0.remove(key.0 as usize);
        }
    }

    /// Add node to tree under `parent`. If `parent` is `None`, the current node is inserted as a root.
    /// It's possible to insert multiple roots into one tree.
    ///
    /// Returns `Err` variant if index is out of the range of `NodeId`.
    pub fn add_node(
        &mut self,
        value: T,
        parent: NodeId,
    ) -> Result<NodeId, std::num::TryFromIntError> {
        Ok(NodeId(u32::try_from(
            self.0.insert(Node { value, parent }),
        )?))
    }

    /// Returns an _inclusive_ iterator over the data associated with the ancestors of a given node ID
    pub fn ancestors(&self, id: NodeId) -> Ancestors<T> {
        Ancestors {
            tree: self,
            state: id,
        }
    }

    /// Returns the number of elements in the tree
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns if the tree is empty
    pub fn is_empty(&self) -> bool {
        self.0.len() == 1
    }

    /// Returns the number of elements the Tree can hold without reallocating
    pub fn capacity(&self) -> usize {
        self.0.capacity()
    }

    /// Clears and initializes the tree,
    /// keeps the allocated memory for reuse. It's necessary to call
    /// this method before using the tree to obtain the `NodeId` of the
    /// root node.
    pub fn clear(&mut self) -> NodeId {
        self.0.clear();
        self.add_internal_root();
        NodeId(0)
    }
}

/// The node struct is not exposed at the API surface
struct Node<T> {
    value: T,
    parent: NodeId,
}

/// Inclusive iterator over the data associated with the ancestors of a given node ID
pub struct Ancestors<'a, T> {
    tree: &'a Tree<T>,
    state: NodeId,
}

impl<'a, T> Iterator for Ancestors<'a, T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.state.0 != 0 {
            if let Some(node) = self.tree.0.get(self.state.0 as usize) {
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

        let root_id = tree.clear();
        let first_id = tree.add_node(value, root_id).unwrap();
        let second_id = tree.add_node(value + 1, first_id).unwrap();
        let third_id = tree.add_node(value + 2, second_id).unwrap();
        let fourth_id = tree.add_node(value + 3, third_id).unwrap();

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
        let root_id = tree.clear();
        tree.add_node(value, root_id).unwrap();

        assert_eq!(tree.len(), 2);
        assert_eq!(tree.len(), tree.0.len());

        tree.add_node(value, root_id).unwrap();

        assert_eq!(tree.len(), 3);
        assert_eq!(tree.len(), tree.0.len());
    }

    #[test]
    fn test_capacity() {
        let mut tree: Tree<usize> = Tree::with_capacity(5);
        tree.clear();
        assert_eq!(tree.capacity(), 5);
        assert_eq!(tree.capacity(), tree.0.capacity());
    }

    #[test]
    fn test_clear() {
        let value = 15;
        let mut tree = Tree::new();

        let root_id = tree.clear();
        for _i in 0..1024 {
            tree.add_node(value, root_id).unwrap();
        }

        assert_eq!(tree.len(), 1025);

        let mut parent = root_id;
        for _i in 0..1024 {
            parent = tree.add_node(value, parent).unwrap();
        }

        assert_eq!(tree.len(), 2049);

        let new_root = tree.clear();

        assert_eq!(new_root.0, 0);
        assert_eq!(tree.len(), 1);
    }
}
