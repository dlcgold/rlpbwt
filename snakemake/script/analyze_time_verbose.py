import os


class TimeObj:
    def __init__(self, command=0, version="", user_time=0, sys_time=0, mem=0):
        self.command = command
        self.version = version
        self.user_time = user_time
        self.sys_time = sys_time
        self.mem = mem

    def __repr__(self):
        return f"[{self.command}, {self.version}, " \
               f"{self.user_time}, {self.sys_time}, {self.mem}] "


def main():
    directory = '../results'
    time_obj_list = []
    for filename in os.listdir(directory):
        file = os.path.join(directory, filename)
        if os.path.isfile(file):
            with open(file) as f:
                time_obj = TimeObj()
                lines = f.readlines()
                if lines[0] == "rlpbwt":
                    time_obj.command = "rlpbwt"
                    time_obj.version = lines[1].strip()
                    time_obj.user_time = float(lines[2].strip())
                    time_obj.sys_time = float(lines[3].strip())
                    time_obj.mem = float(lines[4].strip())
                else:
                    time_obj.command = "pbwt"
                    time_obj.version = lines[1].strip()
                    time_obj.user_time = float(lines[2].strip())
                    time_obj.sys_time = float(lines[3].strip())
                    time_obj.mem = float(lines[4].strip())
                time_obj_list.append(time_obj)
    print(time_obj_list)


if __name__ == "__main__":
    main()
