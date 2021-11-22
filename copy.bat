copy input.json \\192.168.31.153\sambashare\home\steve\Downloads\mndo-server
copy mndo-install.sh \\192.168.31.153\sambashare\home\steve\Downloads\mndo-server
copy mndo-j.sh \\192.168.31.153\sambashare\home\steve\Downloads\mndo-server
copy mndo-update.sh \\192.168.31.153\sambashare\home\steve\Downloads\mndo-server
copy build\libs\MNDOParam.jar \\192.168.31.153\sambashare\home\steve\Downloads\mndo-server
copy atom-properties.json \\192.168.31.153\sambashare\home\steve\Downloads\mndo-server

ssh steve@34.66.107.177 -i C:\Users\billi\Documents\stevecyr "sudo ./mndo-update.sh"