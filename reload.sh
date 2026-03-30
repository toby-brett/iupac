
#!/bin/bash

set -e 

git pull origin main

sudo systemctl daemon-reload
sudo systemctl restart dailychemistry.service
sudo nginx -t 
sudo systemctl reload nginx
