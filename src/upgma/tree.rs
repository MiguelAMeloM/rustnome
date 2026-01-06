use crate::sequences::profile::Profile;

pub enum TreeNode {
    Branch {
        dist: f64,
        profile: Vec<Profile>,
        left: Box<TreeNode>,
        right: Box<TreeNode>,
    },
    Leaf(Vec<Profile>, String),
}

impl TreeNode {
    pub fn merge(self, other: TreeNode, dist: f64, new_prof: Vec<Profile>) -> TreeNode {
        TreeNode::Branch {
            profile: new_prof,
            left: Box::new(self),
            right: Box::new(other),
            dist,
        }
    }

    pub fn prof(&self) -> &Vec<Profile> {
        match self {
            TreeNode::Branch { profile, .. } => profile,
            TreeNode::Leaf(prof, _) => prof,
        }
    }
}
