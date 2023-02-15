import os

try:
    # Try to load a dotenv file
    from dotenv import load_dotenv

    load_dotenv()
except ImportError:
    pass


# Slack token for sending
SLACK_TOKEN = os.getenv("SLACK_TOKEN")
