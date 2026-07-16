function idx = find_line_index(from_bus,to_bus)
    branch = table2array(readtable("C:\Users\19712\Desktop\new line data.xlsx",'Sheet','sheet1','Range','A4:B125'));
    idx = find((branch(:,1) == from_bus & branch(:,2) == to_bus) | (branch(:,1) == to_bus & branch(:,2) == from_bus));
end