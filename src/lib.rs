pub mod distributed;
pub mod errors;
pub mod index;
pub mod map;
pub mod utils;

pub const CRATE_NAME: &str = "mapAD";

pub mod build_info {
    include!(concat!(env!("OUT_DIR"), "/built.rs"));

    pub fn get_software_version() -> String {
        let git_string = GIT_COMMIT_HASH.map_or_else(String::new, |git_commit_hash| {
            format!(
                " ({}{})",
                &git_commit_hash[..8],
                if GIT_DIRTY == Some(true) {
                    "-dirty"
                } else {
                    ""
                }
            )
        });

        format!("{PKG_VERSION}{git_string}")
    }
}
