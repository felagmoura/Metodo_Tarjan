#include <iostream>
#include <list>
#include <stack>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>

#define END_INPUT "testes/graph-test-50000.txt"
#define nulo 0

using namespace std;

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
        ifstream arq (END_INPUT);

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
        //p->back().print();
        p->pop_back();
    }

    void incrementaContadorDeBlocos () {
        //cout << "\n";
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

 

int main () {
    //Grafo grafo;
    //grafo.inicializaMetodoTarjam();
    //cout << "\nNUMERO DE COMPONENTES: " << grafo.cont_blocos << "\n" << "TEMPO DE EXECUÇÃO: " << std::chrono::duration_cast<std::chrono::milliseconds>(exeTime).count() << "ms\n";
    return 0;
}