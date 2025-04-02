import os
import shutil
import requests
from bs4 import BeautifulSoup
import argparse
import random

def fetch_and_parse_url(url):
    response = requests.get(url)
    if response.status_code == 200:
        return BeautifulSoup(response.text, 'html.parser')
    else:
        return None

def get_directory_list(url):
    response = requests.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')
    directories = [link.get('href') for link in soup.find_all('a') if link.get('href').endswith('/')]
    return directories

def save_to_file(data, filename):
    with open(filename, 'w') as f:
        for item in data:
            f.write("%s\n" % item)

def download_file(url, directory):
    response = requests.get(url)
    filename = url.split('/')[-1]
    with open(os.path.join(directory, filename), 'wb') as f:
        f.write(response.content)

def recursive_search(url, target_extension, max_depth, current_depth=0):
    if current_depth > max_depth:
        return []

    soup = fetch_and_parse_url(url)
    if soup:
        file_paths = []
        for link in soup.find_all('a'):
            href = link.get('href')
            if href.endswith(target_extension):
                file_paths.append(url + href)
            elif href.endswith('/'):
                subdirectory = url + href
                file_paths.extend(recursive_search(subdirectory, target_extension, max_depth, current_depth + 1))
        return file_paths
    else:
        return []

def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-g', '--keyword', type=str, help='Keyword to filter directories')
    parser.add_argument('-s', '--skip', type=str, help='Keyword to skip')
    args = parser.parse_args()

    url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/"
    if os.path.exists("dir_index.txt"):
        with open("dir_index.txt", "r") as f:
            dir_list = f.read().splitlines()
    else:
        dir_list = get_directory_list(url)
        save_to_file(dir_list, "dir_index.txt")

    if args.keyword:
        filtered_dirs = [directory for directory in dir_list if args.keyword in directory]
        save_to_file(filtered_dirs, "fil_dir_index.txt")
    else:
        filtered_dirs = dir_list

    filtered_dirs = [directory for directory in filtered_dirs if "sp." not in directory]

    if args.skip:
        filtered_dirs = [directory for directory in filtered_dirs if args.skip not in directory]

    if len(filtered_dirs) < 50:
        selected_dirs = filtered_dirs
    else:
        selected_dirs = random.sample(filtered_dirs, 50)

    save_to_file(selected_dirs, "selected_dir.txt")

    os.makedirs("tmp", exist_ok=True)

    for directory in selected_dirs:
        file_urls = recursive_search(url + directory, 'cds_from_genomic.fna.gz', 3)
        if file_urls:
            download_file(file_urls[0], 'tmp')

if __name__ == "__main__":
    main()

