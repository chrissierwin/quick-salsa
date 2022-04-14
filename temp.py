file = open("temp.txt", "w")
file.write("hello")
file.close()


file = open("temp.txt", "r")
print(file.read())
file.close()