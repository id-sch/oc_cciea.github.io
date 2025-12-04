import os


# --create output directory
dir_out = 'dir_test_gd'

# --check if directory exist, if it doesn't then create
try:
    os.makedirs(dir_out)
except OSError:
    if not os.path.isdir(dir_out):
        raise

with open('{}/readme.txt'.format(dir_out), 'w') as f:
    f.write('Create a new text file!')

dir_list_end = os.listdir()
print("END -------------------------------")
print("Files and directories in  :")
print(dir_list_end)
