pub mod distributed;
pub mod errors;
pub mod index;
pub mod map;
pub mod utils;

pub const CRATE_NAME: &str = "mapAD";

pub mod build_info {
    include!(concat!(env!("OUT_DIR"), "/built.rs"));

    pub fn get_software_version() -> String {
        let git_string = if let Some(git_commit_hash) = GIT_COMMIT_HASH {
            format!(
                " ({}{})",
                &git_commit_hash[..8],
                if let Some(true) = GIT_DIRTY {
                    "-dirty"
                } else {
                    ""
                }
            )
        } else {
            String::new()
        };

        format!("{PKG_VERSION}{git_string}")
    }
}
