set -ex

cd ../../lavaan_fitmodel
tar -zcvf lavaan.tar.gz lavaan
cd ../TRAINING/redux10point2

zone=us-central1-f
prefix=ccarey
script_name=cfa2

name=$prefix-$script_name-training
gcloud compute instances create --zone $zone --custom-cpu 20 --custom-memory 120 --image cfa-benchmarking $name &

wait

sleep 5 # give the instances a moment to fully come alive

gcloud compute scp --zone $zone $script_name.r ../../lavaan_fitmodel/lavaan.tar.gz ../../imp1.Rds ../train_withdraw_outliers.txt ../../full.txt muthen1984.Rds $name:. 
