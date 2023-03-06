# Run on AWS
1. Create AWS account.
2. Create an EC2 instance (follow the instructions)
3. Click in AWS on **Connect to instance** and follow the instructions (getting key)
4. Open terminal/Powershell
5. Open Powershell
    ```console
    ssh -i "inbal.pem" ubuntu@ec2-18-221-210-214.us-east-2.compute.amazonaws.com
    ```

5. install pip:
    ```console
    sudo apt update
    sudo apt install python3-pip
    ```

7. install the project
    ```console
    git clone https://github.com/inbalpreuss/DnaStorage
    cd DnaStorage/
    pip3 install -r requirements.txt
    ```

8. create input file
    ```console
    mkdir -p data/testing
    touch data/testing/input_text.dna
    echo "inbal preuss" > data/testing/input_text.dna
    ```

9. run dna storage
    ```console
    python3 -m dna_storage.main
    cat data/testing/simulation_data.9.text_results_file.dna
    ```

#Run Toky Server

1. Open putty
2. Host Name: toky@194.153.101.31 Port:22
3. Password: toky1234
4. Open a directory named "inbal"
    ```console
   mkdir inbal 
   ```
5. clone repository to inbal directory
    ```conslole
    git clone https://github.com/InbalPreuss/DnaStorage.git 
    ``` 

 7. Run program test_dna.py:
    
    From directory DnaStorage:
    1. Run main:
        ```console
        python3 -m dna_storage.main
        ```
    2. Run test_dna:
    
        1. Because the command:
        ```console
        python3 -m tests.test_dna
        ```
        Doesn't work from some reason on toky server, to run test_dna do the following:
        
        We need to copy test_dna.py file to dna_storage folder
        
        From inbal directory:
        ```console
        cp DnaStorage/tests/test_dna.py DnaStorage/dna_storage
        ```
       Then inside DnaStorage directory run:
       ```console
       python3 -m dna_storage.test_dna
       ```
    3. Run distributed:
        ```console
        nohup python3 -m tests.distributed &
       ```
    4. Run plots:
        ```console
        nohup python3 -m dna_storage.plots &
        ```
       
7. nohup is a command that lets your program run while you are disconnected from putty
    ```console
    nohup python3 -m dna_storage.test_dna &
    ```
   
   Save PID: 
   After writing the commend above, the PID of the process will appear, copy the PID so you can control the process later (to stop the process for example)
   
   Kill PID:
   ```console
   kill 14404
   ```
   
 8. To verify that the program is still running after you disconnect from putty
    ```console
    htop
    ```
    For search press F4
    ```console
    python3
    ```
    
    If you see the program you run you succeeded!! 
    Else, try agin
    
9. Move files from a distant server:
 
    download Winscp 
   

   
