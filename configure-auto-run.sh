echo "making directories......"
folder_path=$(python3 -m src.data_extraction.ConfigGenerator_Deterministic "default_configuration.json" "contour_data/double_variant/")

echo "Folder path: $folder_path"

PYTHON_SCRIPT="src.data_extraction.DataCollector_Deterministic"

run_cases(){
    for ((i=0; i<20; i++)); do
        echo "$folder_path/Relationship_Weight_$i" > "auto_run_files/load_directory.$i.txt"
    	screen -S "$i.contourscreen" -dm
    	
    	screen -S "$i.contourscreen" -X stuff "longjob -28day -c ./parse_screen.sh\n"
    done
}

run_cases

echo "end"

