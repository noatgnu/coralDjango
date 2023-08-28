import re

def replace_special_with_dot(s: str):
    return re.sub(r'[^a-zA-Z0-9.]', '.', s)