copy molecules.txt \\192.168.31.153\sambashare\home\steve\Downloads\mndo-server
copy param-numbers.csv \\192.168.31.153\sambashare\home\steve\Downloads\mndo-server
copy params.csv \\192.168.31.153\sambashare\home\steve\Downloads\mndo-server
copy mndo-install.sh \\192.168.31.153\sambashare\home\steve\Downloads\mndo-server
copy mndo-j.sh \\192.168.31.153\sambashare\home\steve\Downloads\mndo-server
copy mndo-update.sh \\192.168.31.153\sambashare\home\steve\Downloads\mndo-server
copy build\libs\MNDOParam.jar \\192.168.31.153\sambashare\home\steve\Downloads\mndo-server\MNDOParam.jar

ssh steve@34.67.141.141 -i C:\Users\billi\Documents\stevecyr "sudo ./mndo-update.sh"