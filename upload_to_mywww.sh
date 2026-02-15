#!/bin/bash
# Upload the figures to mywww (my own Caesar web server), using scp
# Password authentication is not set up in any special way, so you will be prompted for the password of the mywww user on the server.

# Upload the figs directory under a different name, asked as a command line argument (e.g., "figs_2024_06_01")
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <remote_dir_name>"
    exit 1
fi
REMOTE_DIR="$1"

# Define the local directory containing the figures
LOCAL_DIR="figs"

# Define the remote directory path on the server
REMOTE_PATH="public_html/EPOS_UrQMD/$REMOTE_DIR"
# Use scp to upload the local directory to the remote server
scp -r "$LOCAL_DIR" molnarmatyas@caesar.elte.hu:"$REMOTE_PATH"
# Ready I guess, the figures should now be accessible at http://caesar.elte.hu/~molnarmatyas/$REMOTE_DIR/
