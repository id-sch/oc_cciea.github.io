# hello_world.py
print("Hello, World!  Hello, World!")

file = open("hello_world.org", "w")
file.write("* Hello World\n+ this is my first file created using Python!\n+ you need to set read/write permissions in Actions/General on Github\n+ I don't know why this is not working now.")
file.close()
print("File created successfully!")
