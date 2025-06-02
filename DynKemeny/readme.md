## Environment
- Ubuntu 22.04
- GCC 11.2.0

## Compile
- make

## Run
Run an algorithm to calculate Kemeny constant.
```
main -d <dataset_path> -a <algorithm_name> -static [options]
```
Descriptions:
- algorithm_name: ttf, spantree, forestmc, lewalk
- options:
    - eps: relevant to the index size. (0.5 by default)
    - static: running on static algorithms (read graph file "graph.txt").
    - dynamic: running on dyanmic algorithms (read graph file "base_graph.txt" and updates file "add_edges.txt", "del_edges.txt")
        - imp: using ImprovedSM algorithm (have to use "-a ttf")
        - basic: using BasicSM algorithm (have to use "-a ttf")
    - log: the log file path.
Example:
```sh
# calculate kemeny constant for static grpah "powergrid" with epsilon = 0.1 using TTF
./main -d powergrid -a ttf -eps 0.1 -static
# calculate kemeny constant for dynamic grpah "powergrid" with 100 edges deleted and 100 edges inserted using ImprovedSM
./main -d powergrid -a ttf -dynamic -imp -log result_dynamic.txt
```

## Run experiment
```sh
sh shell/exp-static.sh
```