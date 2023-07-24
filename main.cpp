#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <stack>
#include <set>
#include <algorithm>
#include <chrono>
#include <list>

#define FILE_PATH "testes/grafo_10V_20E_3A.txt"
#define nulo 0

using namespace std;

class SimpleGraph
{
public:
  [[maybe_unused]] explicit SimpleGraph(int vertex_num)
  {
    this->vertices = std::vector<Vertex>(vertex_num + 1);
    for (int i = 1; i < vertex_num; i++)
    {
      vertices[i] = Vertex();
    }
  }

  explicit SimpleGraph(const std::string &file_path)
  {
    std::ifstream input_file(file_path);
    ParseVerticesFromFile(input_file);
    ParseEdgesFromFile(input_file);

    input_file.close();
  }

#pragma region main_public_methods

  std::vector<std::set<int>> FindBridges(){
    std::vector<std::set<int>> bridges = std::vector<std::set<int>>(vertices.size());
    std::cout << "Brigdes: ";
    for (int i = 1; i < vertices.size(); i++)
    {
      for (int j = 0; j < vertices[i].neighbors.size(); j++)
      {
        if(isBridge(i, vertices[i].neighbors[j])){
          std::set<int> bridge = std::set<int>();
          bridge.insert(i);
          bridge.insert(vertices[i].neighbors[j]);
          bridges.push_back(bridge);
          //std::cout << "(" << i << ", " << vertices[i].neighbors[j] << ") ";

        }
      }
      
    }
    return bridges;
  }

  std::vector<std::vector<int>> getBridgesFromVertex(int vertex_origin){
    std::vector<std::vector<int>> bridges = std::vector<std::vector<int>>();

    for (int i = 0; i < vertices[vertex_origin].neighbors.size(); i++)
    {
      //std::cout << "testando para " << vertex_origin << " e " << vertices[vertex_origin].neighbors[0] << "\n";
      //Pega o primeiro item da lista de vizinhos do vertice para testar se é ponte
      if(isBridge(vertex_origin, vertices[vertex_origin].neighbors[0])){
        //std::cout << "EH PONTEEE\n";
        //Caso seja ponte, adiciona na lista de pontes
        std::vector<int> bridge = std::vector<int>(2);
        bridge[0] = vertex_origin;
        bridge[1] = vertices[vertex_origin].neighbors.back();   //Para testar se é ponte, removemos e então reinserimos o item testado na lista de vizinhos, no final
        bridges.push_back(bridge);
      }

      
    }
    return bridges;
  }

  void dfs()
  {
    std::vector<bool> visited = std::vector<bool>(vertices.size(), false);
    for (int i = 1; i < vertices.size(); i++)
    {
      if (!visited[i])
      {
        dfs(i, visited);
      }
    }
  }

  std::vector<std::set<int>> GetBlocksFromDisjointPaths()
  {
    int blocks_count = 0;
    std::vector<std::set<int>> blocks = std::vector<std::set<int>>();

    for (int i = 1; i < vertices.size(); i++)
    {
      blocks.emplace_back();
      for (int j = 1; j < vertices.size(); j++)
      {
        std::cout << "checking for " << i << " and " << j << "\n";
        if (CheckDifferentPaths(DfsForCheckingAllPaths(i, j)))
        {
          blocks[blocks_count].insert(i);
          blocks[blocks_count].insert(j);
        }
      }
      if (blocks[blocks_count].empty())
      {
        for (auto adjacent_vertex : vertices[i].neighbors)
        {
          blocks[blocks_count].insert(i);
          blocks[blocks_count].insert(adjacent_vertex);
          blocks_count++;
          blocks.emplace_back();
        }
      }
      else
      {
        blocks_count++;
      };
    }

    std::set<std::set<int>> unique_paths(blocks.begin(), blocks.end());

    DeleteRedundantPaths(unique_paths);
    CheckForBlocksOfTwo(unique_paths);
    PrintBlocks(unique_paths);

    return blocks;
  }

  static void PrintBlocks(const std::set<std::set<int>> &paths)
  {
    std::cout << std::endl;
    for (const auto &path : paths)
    {
      for (int vertex : path)
      {
        std::cout << vertex << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "Temos um total de " << paths.size() << " blocos." << std::endl;
  }

  void printGraph()
  {
    for (int i = 1; i < vertices.size(); i++)
    {
      std::cout << i << ": ";
      for (int j = 0; j < vertices[i].neighbors.size(); j++)
      {
        std::cout << vertices[i].neighbors[j] << " ";
      }
      std::cout << std::endl;
    }
  }


  std::vector<std::set<int>> GetBlocksFromArticulations()
  {
    std::vector<bool> articulations = FindArticulations();

    //imprime as articulações
    //std::cout << "Articulacoes: ";
    for (int i = 1; i < articulations.size(); i++)
    {
      if (articulations[i])
      {
        //std::cout << i << " ";
      }
    }

    std::vector<std::set<int>> blocks = std::vector<std::set<int>>();

    //ordena os vertices para priorizar os vizinhos não articulações
    std::vector<int> sorted_vertices = std::vector<int>();
    for (int i = 1; i < vertices.size(); i++)
    {
      //Recupera todas as nao articulações do vertice
      std::vector<int> nonArt = std::vector<int>();
      for (int &neighbor : vertices[i].neighbors)
      {
        if (!articulations[neighbor])
        {
          nonArt.push_back(neighbor);
        }
      }
      for (int &neighbor : vertices[i].neighbors)
      {
        if (articulations[neighbor])
        {
          nonArt.push_back(neighbor);
        }
      }
      vertices[i].neighbors = nonArt;
    }
    
    //Faz a busca de blocos para cada articulação
    for (int i = 1; i < vertices.size(); i++)
    {
      if (articulations[i])
      {
        std::vector<std::set<int>> blocksFromArt = GetBlocksByArticulation(i, articulations);
        for (std::set<int> &block : blocksFromArt)
        {
          blocks.push_back(block);
        }
      }
    }

    //Remove blocos repetidos
    std::set<std::set<int>> unique_blocks(blocks.begin(), blocks.end());
    blocks = std::vector<std::set<int>>(unique_blocks.begin(), unique_blocks.end());

    //blocks = GetBlocksByArticulation(v, articulations);
    return blocks;
  }

  std::vector<std::set<int>> GetBlocksByArticulation(int root, std::vector<bool> articulations)
  {
    std::vector<int> call_stack = std::vector<int>();
    std::vector<int> uniqueBlock = std::vector<int>();
    std::vector<std::set<int>> blocks = std::vector<std::set<int>>();
    //Vetores para a busca em profundidade
    std::vector<int> discoverTime = std::vector<int>(vertices.size(), 0);
    std::vector<int> father = std::vector<int>(vertices.size(), 0);
    //std::cout << "Root: " << root << "\n";

    int counter = 1;    //Contador para o tempo de descoberta
    int rootNeighbors = 0;  //Contador para o numero de vizinhos da raiz, para garantir que visitará todos
    discoverTime[root] = counter;
    counter++;

    bool found_unvisited_neighbor = false;
    int current_vertex;

    std::vector<std::vector<int>> bridgesOfVertex = getBridgesFromVertex(root);   //Encontra todas as pontes conectadas à raiz
    //Insere todas as pontes do vértice raiz, pois são blocos
    for (int i = 0; i < bridgesOfVertex.size(); i++)
    {
      blocks.push_back(std::set<int>());
      for (int j = 0; j < bridgesOfVertex[i].size(); j++)
      {
        blocks[blocks.size() - 1].insert(root);
        blocks[blocks.size() - 1].insert(bridgesOfVertex[i][j]);
        rootNeighbors++;    //Incrementa o numero de vizinhos da raiz por cada ponte
        //Marca o tempo de descoberta como -1 para ignorar na busca em profundidade
        discoverTime[bridgesOfVertex[i][j]] = -1;
      }
    }

    //Insere a raiz na pilha de chamadas e no bloco unico
    call_stack.push_back(root);
    uniqueBlock.push_back(root);
    
    while (call_stack.size() > 1 || rootNeighbors < vertices[root].neighbors.size()){
      current_vertex = call_stack.back();
      //std::cout << "Current: " << current_vertex << "\n";
      
      if(current_vertex <= 0 || current_vertex >= vertices.size())
      {
        //std::cout << "Erro: current_vertex < 0\n";
        break;
      }
      for (int i = 0; i < vertices[current_vertex].neighbors.size(); i++)
      {
        
        int neighbor = vertices[current_vertex].neighbors[i];
        if(neighbor <= 0 || neighbor >= vertices.size())
        {
          //std::cout << "Erro: neighbor " << neighbor << "\n";
          break;
        }
        //Ignora os vizinhos ja descobertos anteriormente
        else if (discoverTime[neighbor] < 0){
          continue;
        }
        // Aresta de árvore
        else if (discoverTime[neighbor] == 0)
        {
          //std::cout << "Arvore: " << current_vertex << " -> " << neighbor << "\n";
          call_stack.push_back(neighbor);   //Adiciona o vizinho na pilha de chamadas
          uniqueBlock.push_back(neighbor);  //Adiciona a aresta ao bloco
          father[neighbor] = current_vertex;  //Define o pai do vizinho como o atual
          discoverTime[neighbor] = counter;  //Define o tempo de descoberta do vizinho
          counter++;
          found_unvisited_neighbor = true;

          if(call_stack.back() < 0){
              //std::cout << "Erro: call_stack inseriu " << call_stack.back() <<"\n";
              while(true);
              break;
            }
          break;
        }
        // aresta de retorno
        else if (discoverTime[neighbor] < discoverTime[current_vertex] && neighbor != father[current_vertex])
        {
          //std::cout << "Retorno: " << current_vertex << " -> " << neighbor << "\n";
          //Se colidiu com uma aresta de articulação que não foi a raiz, pertence a outro bloco
          if (articulations[neighbor] && neighbor != root)
          {
            //std::cout << "Colidiu com: " << current_vertex << " retornando para " << neighbor << "\n";
            //Remove os vértices até achar a articulação que colidiu
            while (call_stack.back() != neighbor && call_stack.size() > 0)
            {
              //std::cout << "Removendo " << call_stack.back() << " procurando: " << neighbor << "\n";
              if(call_stack.back() < 0){
                //std::cout << "Erro: call_stack.back() " << call_stack.back() << " size stack: " << call_stack.size() << "\n";
                while(true);
                break;
              }
              if(call_stack.back() == root){
                break;
              }
              call_stack.pop_back();
              uniqueBlock.pop_back();
            }
            current_vertex = call_stack.back();
            //std::cout << "Compara " << call_stack.back() << " com " << neighbor << "\n";
          }
          
        }

        if(call_stack.size() == 0){
          //std::cout << "Erro: call_stack.size() == 0\n";
          break;
        }

      }

      if(call_stack.size() == 0){
          //std::cout << "Erro: call_stack.size() == 0\n";
          break;
        }

      //std::cout << "debug\n";
      if (!found_unvisited_neighbor)
      {
        //std::cout << "Backtrack: " << current_vertex << "\n";
        call_stack.pop_back();  //Remove o vértice de "ponta solta" da pilha de chamadas, pois é uma ponte que já foi tratada
        if(father[current_vertex] == call_stack.back() && vertices[current_vertex].neighbors.size() == 1){
          uniqueBlock.pop_back();
         
        }
        //Se retornou a raiz, encontrou um bloco
        if (call_stack.back() == root)
        {
          //Insere o bloco na lista de blocos
          blocks.push_back(std::set<int>());

          //Zera os tempos de descoberta
          for (int i = 0; i < vertices.size(); i++)
          {
            if(discoverTime[i] > 0)
              discoverTime[i] = 0;
          }
          //Negativa os tempos de descoberta dos vertices ja encontrados anteriormente, pertencentes ao bloco
          for(int & v : uniqueBlock){
            if(IsAdjacent(root, v)){
              rootNeighbors++;
            }
            discoverTime[v] = -1;
          }
          while (!uniqueBlock.empty())
          {
            blocks[blocks.size() - 1].insert(uniqueBlock.back());
            uniqueBlock.pop_back();
          }
          //std::cout << "Inseriu o bloco\n";
          uniqueBlock.push_back(root);
        }
      }
      found_unvisited_neighbor = false;
      //std:: cout << "root neighbors: " << rootNeighbors << " real: " << vertices[root].neighbors.size() << " s: " << call_stack.size() <<"\n";
    }

    return blocks;
  }

#pragma endregion

private:
  struct Vertex
  {
    std::vector<int> neighbors;

    explicit Vertex(std::vector<int> neighbors = std::vector<int>()) : neighbors(std::move(neighbors)) {}

    void AddNeighbor(int successor_id)
    {
      this->neighbors.push_back(successor_id);
    }
  };
  std::vector<Vertex> vertices;

#pragma region parse_file_to_graph
  void ParseEdgesFromFile(std::istream &input_file)
  {
    std::string line;
    int value1, value2;
    while (std::getline(input_file, line))
    {
      std::stringstream ss(line);
      if (ss >> value1 >> value2)
      {
        this->AddEdge(value1, value2);
      }
      else
      {
        std::cerr << "Linha com formato invalido: " << line << std::endl;
      }
    }
  }

  void ParseVerticesFromFile(std::istream &input_file)
  {
    std::string line;
    std::getline(input_file, line);
    std::stringstream ss(line);
    int value1, value2;
    if (ss >> value1 >> value2)
    {
      this->vertices = std::vector<Vertex>(value1 + 1);
      for (int i = 1; i < value1; i++)
      {
        vertices[i] = Vertex();
      }
    }
    else
    {
      std::cerr << "Linha com formato invalido: " << line << std::endl;
    }
  }

  void AddEdge(int origin_id, int destiny_id)
  {
    vertices[origin_id].AddNeighbor(destiny_id);
    vertices[destiny_id].AddNeighbor(origin_id);
  }
#pragma endregion

#pragma region detect_blocks_checking_for_biconnected_vertices
  std::vector<std::vector<int>> DfsForCheckingAllPaths(int start_vertex, int end_vertex)
  {
    std::vector<std::vector<int>> paths;
    std::stack<std::pair<int, std::vector<int>>> stack;
    std::vector<int> initial_path = {start_vertex};
    stack.push({start_vertex, initial_path});

    while (!stack.empty())
    {
      int current_vertex = stack.top().first;
      std::vector<int> current_path = stack.top().second;
      stack.pop();

      if (current_vertex == end_vertex)
      {
        paths.push_back(current_path);
        if (CheckDifferentPaths(paths))
        {
          return paths;
        }
      }
      else
      {
        for (int neighbor : vertices[current_vertex].neighbors)
        {
          bool already_visited = std::find(current_path.begin(), current_path.end(), neighbor) != current_path.end();
          if (!already_visited)
          {
            std::vector<int> new_path = current_path;
            new_path.push_back(neighbor);
            stack.push({neighbor, new_path});
          }
        }
      }
    }
    return paths;
  }

  static bool CheckDifferentPaths(std::vector<std::vector<int>> paths)
  {
    bool flag = true;
    for (int i = 0; i < paths.size() - 1; i++)
    {
      for (int j = i + 1; j < paths.size(); j++)
      {
        std::set<int> set1(paths[i].begin(), paths[i].end());
        for (int element : paths[j])
        {
          if (set1.count(element) > 0 && element != paths[i][0] && element != paths[i][paths[i].size() - 1])
          {
            flag = false;
            break;
          }
        }
        if (flag)
        {
          return flag;
        }
        flag = true;
      }
    }
    return false;
  }

  std::set<std::set<int>> DeleteRedundantPaths(std::set<std::set<int>> &paths)
  {
    for (auto it = paths.begin(); it != paths.end();)
    {
      if (it->empty())
      {
        it = paths.erase(it);
      }
      else
      {
        ++it;
      }
    }
    for (const auto &subset1 : paths)
    {
      for (const auto &subset2 : paths)
      {
        if (subset1 < subset2 && std::includes(subset2.begin(), subset2.end(), subset1.begin(), subset1.end()))
        {
          paths.erase(subset2);
          break;
        }
      }
    }

    auto it = paths.begin();
    while (it != paths.end())
    {
      auto jt = it;
      ++jt;
      bool erase_it = false;
      while (jt != paths.end())
      {
        int count = 0;
        for (auto x : *it)
        {
          if (jt->count(x) > 0)
          {
            count++;
          }
        }
        if (count >= 2)
        {
          if (it->size() > jt->size())
          {
            erase_it = true;
            break;
          }
          else
          {
            jt = paths.erase(jt);
          }
        }
        else
        {
          ++jt;
        }
      }
      if (erase_it)
      {
        it = paths.erase(it);
      }
      else
      {
        ++it;
      }
    }

    return paths;
  }

  std::set<std::set<int>> CheckForBlocksOfTwo(std::set<std::set<int>> &paths)
  {
    std::vector<std::set<int>> new_blocks = std::vector<std::set<int>>();
    int count = 0;
    for (const auto &set1 : paths)
    {
      for (const auto &set2 : paths)
      {
        for (int element : set1)
        {
          for (int element2 : set2)
          {
            if (IsAdjacent(element, element2) && set1.count(element2) == 0 && set2.count(element) == 0)
            {
              new_blocks.emplace_back();
              new_blocks[count].insert(element);
              new_blocks[count].insert(element2);
              count++;
            }
          }
        }
      }
    }
    for (const auto &new_block : new_blocks)
    {
      bool flag = false;
      for (const auto &block : paths)
      {
        if (std::includes(block.begin(), block.end(), new_block.begin(), new_block.end()))
        {
          flag = true;
        }
      }
      if (!flag)
      {
        paths.insert(new_block);
      }
    }
    return paths;
  }

  bool IsAdjacent(int vertex1, int vertex2)
  {
    for (auto vertex : vertices[vertex1].neighbors)
    {
      if (vertex == vertex2)
      {
        return true;
      }
    }
    return false;
  }

#pragma endregion

  void dfs(int vertex_id, std::vector<bool> &visited)
  {
    std::vector<int> call_stack = std::vector<int>();
    bool found_unvisited_neighbor = false;
    int current_vertex;
    call_stack.push_back(vertex_id);

    while (!call_stack.empty())
    {
      current_vertex = call_stack.back();
      visited[current_vertex] = true;

      for (int &neighbor : vertices[current_vertex].neighbors)
      {
        if (!visited[neighbor])
        {
          call_stack.push_back(neighbor);
          found_unvisited_neighbor = true;
          break;
        }
      }
      if (!found_unvisited_neighbor || vertices[current_vertex].neighbors.empty())
      {
        call_stack.pop_back();
      }
      found_unvisited_neighbor = false;
    }
  }

  void findingBlocksByArt(int vertex_id, std::vector<bool> &visited)
  {
    std::stack<int> call_stack = std::stack<int>();
    bool found_unvisited_neighbor = false;
    int current_vertex;
    call_stack.push(vertex_id);

    while (!call_stack.empty())
    {
      current_vertex = call_stack.top();
      visited[current_vertex] = true;

      for (int &neighbor : vertices[current_vertex].neighbors)
      {
        if (!visited[neighbor])
        {
          call_stack.push(neighbor);
          found_unvisited_neighbor = true;
          break;
        }
      }
      if (!found_unvisited_neighbor || vertices[current_vertex].neighbors.empty())
      {
        call_stack.pop();
      }
      found_unvisited_neighbor = false;
    }
  }

  bool isArticulation(int vertex_id)
  {
    std::vector<bool> visited = std::vector<bool>(vertices.size(), false);
    visited[vertex_id] = true;
    int dfsRoot = vertex_id != 1 ? 1 : 2;
    dfs(dfsRoot, visited);
    for (int i = 1; i < vertices.size(); i++)
    {
      if (!visited[i])
      {
        return true;
      }
    }
    return false;
  }

  std::vector<bool> FindArticulations()
  {
    std::vector<bool> articulations = std::vector<bool>(vertices.size(), false);
    for (int i = 1; i < vertices.size(); i++)
    {
      if (isArticulation(i))
      {
        articulations[i] = true;
      }
    }
    return articulations;
  }

  void removeEdge(int v, int w){
    auto destIndex = std::find(this->vertices[v].neighbors.begin(), this->vertices[v].neighbors.end(), w);
    auto originIndex = std::find(this->vertices[w].neighbors.begin(), this->vertices[w].neighbors.end(), v);
    if(destIndex != this->vertices[v].neighbors.end() && originIndex != this->vertices[w].neighbors.end()){
      this->vertices[v].neighbors.erase(destIndex);
      this->vertices[w].neighbors.erase(originIndex);
    }
  }

  bool isBridge(int vertex_origin, int vertex_destination){
    std::vector<bool> visited = std::vector<bool>(vertices.size(), false);
    bool is_bridge = false;
    //Remove a aresta para testar se é ponte
    removeEdge(vertex_origin, vertex_destination);  

    dfs(vertex_origin, visited);  //Realiza a busca, se algum vertice não for visitado, a aresta é ponte
    //Verifica se algum vertice não foi visitado
    for(int i = 1; i < vertices.size(); i++){
      if(!visited[i]){
        is_bridge = true;
        //std::cout << "A aresta " << vertex_origin << " - " << vertex_destination << " e ponte" << std::endl;
      }
    }

    //Retorna a aresta removida para o grafo, mas insere no final do vetor de vizinhos de cada vértice
    AddEdge(vertex_origin, vertex_destination);

    return is_bridge;
  }

  
};

/******************************************************************************/
/*                                  ARESTA                                    */
/******************************************************************************/

class Aresta {
 public:
    int v;
    int w;

    Aresta (int v, int w) {
        this -> v = v;
        this->w = w;
    }

    void print () {
        cout << "(" << v << ", " << w << "); ";
    }
};

/******************************************************************************/
/*                                  GRAFO                                     */
/******************************************************************************/

class Grafo {
    list<int>* list_adj;
    int num_vertices;
    int num_arestas;
    
 public:
    int cont_blocos;

    /******************************************************************************/
    /* CONSTRUTOR DO GRAFO                                                        */
    /*******************************************************************************
     * - v -> vértice visitado;
     * - w -> adjacente do vértice;
     * - arq -> arquivo de entrada de dados do grafo;
     * _____________________________________________________________________________
     * 
     * inicializa um objeto grafo com os dados lidos de um arquivo .txt;
     * iniicializa e preenche uma lista adjacente, que serve como representação 
     * do grafo;
     * 
    *******************************************************************************/
    Grafo (int v = 0, int w = 0) {
        ifstream arq (FILE_PATH);

        if (!arq) 
            throw runtime_error ("arquivo não encontrado");
        
        cont_blocos = 0;
        arq >> num_vertices >> num_arestas;
        list_adj = new list<int>[++num_vertices];

        while (arq >> v && arq >> w) {
            list_adj[v].push_back(w);
            list_adj[w].push_back(v);
        }

        arq.close();
    }

    /******************************************************************************/
    /* IMPRIME A LISTA ADJACENTE                                                  */
    /*******************************************************************************
     * - v -> vértice cabeça da lista adjacente;
     * - it -> percorre pelos adjacentes de v;
     * _____________________________________________________________________________
     * 
     * v percorre por cada vértice do grafo;
     * enquanto it imprime os vértices adjacentes a v; 
     * 
    *******************************************************************************/
    void imprimeLista () {
        list<int>::iterator it;
        
        for (int v = 0; v < num_vertices; v++) {
            cout << "["<< v <<"] -> ";
            for (it = list_adj[v].begin(); it != list_adj[v].end(); ++it)
                cout << *it << "; "; 
            cout << "\n";
        }
    }

    void imprimeRemoveDaPilha (list<Aresta>* p) {
        p->back().print();
        p->pop_back();
    }

    void incrementaContadorDeBlocos () {
        cout << "\n";
        cont_blocos++;
    }

    /******************************************************************************/
    /* VERIFICA SE O VERTICE É UMA ARTICULAÇÃO                                    */
    /*******************************************************************************
     * - tv -> tempo de descoberta do vértice;
     * - lw -> low value do vértice adjacente;
     * - f -> numero de filhos do vértice
     * _____________________________________________________________________________
     * 
     * true quando o vértice for uma articulação;
     * 
    *******************************************************************************/
    bool articulacaoEncontrada (int tv, int lw, int f) {
        return (tv == 1 && f > 1) || (tv > 1 && lw >= tv);
    }

    /******************************************************************************/
    /* VERIFICA SE UM BLOCO FOI ENCONTRADO                                        */
    /*******************************************************************************
     * - p -> pilha dos vértices visitados;
     * - v -> vértice de articulação;
     * - w -> adjacente do vértice;
     * _____________________________________________________________________________
     * 
     * true quando a articulação é encontrada na pilha;
     * neste momento, todas as arestas que estiverem na pilha formam um 
     * componete biconexo;
     * 
    *******************************************************************************/
    bool blocoEncontrado (list<Aresta>* p, int v, int w) {
        return (p->back().v != v || p->back().w != w);
    }

    /******************************************************************************/
    /* INICIALIZA E EXECUTA A DFS DE TARJAM                                       */
    /*******************************************************************************
     * - td[] -> tempo de descoberta dos vértices;
     * - lv[] -> low value, o valor do td mínimo de um vértice alcançável;
     * - pai[] -> o pai do vértice visitado;
     * - pilha[] -> pilha de vértices que fazem parte do componete biconectado;
     * _____________________________________________________________________________
     * 
     * inicializa as variáveis;
     * executa a busca enquanto existir vértices que ainda não foram descobertos;
     * esvazia a pilha se ela não estiver vazia;
     * incrementa o contador de blocos se a pilha não estava vazia depois da busca;
     * 
    *******************************************************************************/
    void inicializaMetodoTarjam () {
        static int* td = new int[num_vertices];
        int* lv = new int[num_vertices];
        int* pai = new int[num_vertices];

        list<Aresta>* pilha = new list<Aresta>[num_arestas];
        bool pilha_vazia;

        for (int i = 0; i < num_vertices; i++) {
            td[i] = nulo;
            lv[i] = nulo;
            pai[i] = nulo;
        }

        for (int v = 0; v < num_vertices; v++) {
            if (td[v] == nulo)
                dfsTarjam (v, td, lv, pilha, pai);

            pilha_vazia = false;

            while (pilha->size() > 0) {
                imprimeRemoveDaPilha(pilha);
                pilha_vazia = true;
            }

            if (pilha_vazia)
                incrementaContadorDeBlocos();
        }
    }

 private:

    /******************************************************************************/
    /* FUNÇÃO RECURSIVE DE DFS ADAPTADA POR TARJAM                                       */
    /*******************************************************************************
     * - t -> tempo global;
     * - filho -> contador de filhos de um vértice;
     * _____________________________________________________________________________
     * 
     * inicializa as variáveis de tempo de descoberta e low value 
     * do vertice visitado;
     * para todo vértice adjacente, verificar se ele já foi descoberto: 
     *  SE NÂO,
     *      insere essa aresta na pilha,
     *      executa a busca com o adjacente como argumento,
     *      no callback da função atualizaa os low values de v e w,
     *      procura pela articulação na pilha para encontrar os 
     *      componentes biconexos,
     *  SE SIM,
     *      atulizar os low values de v e w;
     *      inserir essa aresta na pilha;
     *      
    *******************************************************************************/
    void dfsTarjam (int v, int td[], int lv[], list<Aresta>* pilha, int pai[]) {
        static int t = 0;
        int filhos = 0;

        td[v] = lv[v] = ++t;
        
        list<int>::iterator it;
        for (it = list_adj[v].begin(); it != list_adj[v].end(); it++) {
            int w = *it;

            if (td[w] == nulo) {
                filhos++;
                pai[w] = v;
                pilha->push_back(Aresta(v, w));
                
                dfsTarjam (w, td, lv, pilha, pai);
                lv[v] = min(lv[v], lv[w]);
                
                if (articulacaoEncontrada (td[v], lv[w], filhos)) {
                    cout << "Arestas do Bloco " << cont_blocos + 1 << ": ";
                    while (blocoEncontrado (pilha, v, w)) {
                        imprimeRemoveDaPilha(pilha);
                    }
                    imprimeRemoveDaPilha(pilha);
                    incrementaContadorDeBlocos();
                }

            } else if (w != pai[v]) {
                lv[v] = min(lv[v], td[w]);
                
                if (td[w] < td[v]) 
                    pilha->push_back(Aresta(v, w));
            }
        }
    }
};

#pragma region input_methods
std::string InputFilePath()
{
  std::string input;
  std::cout << "\nInsira o caminho do arquivo de grafos: ";
  std::cin >> input;
  return input;
}

#pragma endregion

void geraHUD () {
  cout 
    << "************************************************************************************\n"
    << "Selecione o metodo de identificacao de blocos: \n"
    << "[1] buscar por caminhos disjuntos entre pares de vertices\n"
    << "[2] buscar por articulacoes\n"
    << "[3] metodo de Tarjam\n"
    << "[0] encerrar programa\n" 
    << "************************************************************************************\n";
}

void metodoCaminhosDisjuntos () {
  SimpleGraph g(FILE_PATH);
  cout << "Iniciando...\n";
  std::vector<std::set<int>> block = g.GetBlocksFromDisjointPaths ();
  std::cout << "Tamanho do bloco: " << block.size() << std::endl;
}

void metodoIdentificacaoPorArticulacao () {
  SimpleGraph g(FILE_PATH);
  cout << "Iniciando...\n";
  std::vector<std::set<int>> block = g.GetBlocksFromArticulations ();

  for (int i = 0; i < block.size(); i++) {
    std::cout << "Bloco: " << i + 1 << ": ";
    for (int j : block[i]) {
      std :: cout << j << " ";
    }
    std::cout << std::endl;
  }

  std::cout << "Tamanho do bloco: " << block.size() << std::endl;
}

void metodoTarjam () {
  Grafo grafo;
  grafo.inicializaMetodoTarjam ();
  cout << "NUMERO DE BLOCOS: " << grafo.cont_blocos << "\n";
}

bool executarMetodoSelecionado (bool fim = false) {
  int op;
  cout << "Metodo Selecionado: ";
  cin >> op;
  cout << "************************************************************************************\n";

  switch (op)
  {
  case 1:
    metodoCaminhosDisjuntos ();
    break;
  
  case 2:
    metodoIdentificacaoPorArticulacao ();
    break;

  case 3:
    metodoTarjam ();
    break;
  
  case 0:
  default:
    fim = true;
    cout 
      << "PROGRAMA ENCERRADO\n"
      << "************************************************************************************\n";
    break;
  }

  return fim;
}


int main()
{
  bool fim = false;
  while (!fim) {
    geraHUD();
    fim = executarMetodoSelecionado();
  }
  
  return 0;
}


