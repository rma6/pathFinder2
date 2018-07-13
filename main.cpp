#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <exception>
#include <chrono>
#include "omp.h"

#define DEFAULT_MAX_THREADS 4

using namespace std;

class coords
{
    public:
    int line, column, direction;

    vector<vector<int>> reg;
    vector<int> path;
    int cells;

    coords(int a, int b, int c, int d, vector<int>& e, vector<vector<int>>& f) //for walker
    {
        line=a;
        column=b;
        direction=c;
        cells=d;
        path=e;
        reg=f;
    }

    coords(int a=0, int b=0, int l=0, int c=0, int sl=0, int sc=0, int cc=0, vector<int> pre=vector<int>()) //for main
    {
        line=a;
        column=b;
        direction=-1;
        cells=cc;
        path=pre;
        reg=vector<vector<int>>(l, vector<int>(c));
        if(l!=0 && c!=0)
        {
            reg[sl][sc]=-1; //thats for free cell preference
        }
    }

    coords& operator=(const coords& other)
    {
        this->line=other.line;
        this->column=other.column;
        this->direction=other.direction;
        this->reg=other.reg;
        this->path=other.path;
        this->cells=other.cells;
        return *this;
    }
};


void dead_end_thread(coords);
int dead_end_sc(coords ic);
void walker(coords);
void progress();


vector<string> CM;
vector<vector<int>> ACM;
vector<vector<int>> dead_end_vector_position;
vector<vector<int>> dead_ends_paths;
vector<vector<int>> s_reg;
coords start;
vector<int> solution, de_append, de_preppend;//de_append: path to append in case the starting cell is in a dead end
int lines=0, columns=0, open_cells=0, max_threads, runing_threads=0, sc, scide=0; //scide: starting cell in dead_end
bool isrunning=true, dece=false; //dece: dead_end conditional enable
condition_variable_any pcv;

mutex** nav_mtxs;//fixed
mutex m_dep;

int main(int argc, char* argv[])
{
    max_threads=argc>1 ? stoi(argv[1]) : DEFAULT_MAX_THREADS;
    vector<coords> dead_ends; //dead_ends

    //read file
    fstream file;
    file.open("matrix.txt", fstream::in);
    file >> start.line >> start.column;
    for(string temp; getline(file, temp); lines++)
    {
        CM.push_back(temp);
    }
    lines--;//for compensation
    CM.erase(CM.begin());


    //conversion and inicialization
    ACM.resize(lines);
    s_reg.resize(lines);
    dead_end_vector_position.resize(lines);
    nav_mtxs  = (mutex**)malloc(sizeof(mutex*)*lines);
    columns=CM[0].size();
    for(int i=0; i<lines; i++)
    {
        dead_end_vector_position[i].resize(columns);
        nav_mtxs [i] = (mutex*)malloc(sizeof(mutex)*columns);
        ACM[i].resize(columns);
        s_reg[i].resize(columns);
    }
    
    //corners
    if((ACM[0][0] = CM[0][0]=='0' ? (CM[1][0]=='0'? 1:0)+(CM[0][1]=='0'? 1:0) : 0)==1)
    {
        dead_ends.emplace_back(0, 0);
    }
    if(ACM[0][0]!=0)
    {
        open_cells++;
    }
    if((ACM[lines-1][0] = CM[lines-1][0]=='0' ? (CM[lines-2][0]=='0'? 1:0)+(CM[lines-1][1]=='0'? 1:0) : 0)==1)
    {
        dead_ends.emplace_back(lines-1, 0);
    }
    if(ACM[lines-1][0]!=0)
    {
        open_cells++;
    }
    if((ACM[0][columns-1] = CM[0][columns-1]=='0' ? (CM[1][columns-1]=='0'? 1:0)+(CM[0][columns-2]=='0'? 1:0) : 0)==1)
    {
        dead_ends.emplace_back(0, columns-1);
    }
    if(ACM[0][columns-1]!=0)
    {
        open_cells++;
    }
    if((ACM[lines-1][columns-1] = CM[lines-1][columns-1]=='0' ? (CM[lines-2][columns-1]=='0'? 1:0)+(CM[lines-1][columns-2]=='0'? 1:0) : 0)==1)
    {
        dead_ends.emplace_back(lines-1, columns-1);
    }
    if(ACM[lines-1][columns-1]!=0)
    {
        open_cells++;
    }

    //edges        
    for(int i=0, j=1; j<columns-1; j++)
    {
        if((ACM[i][j] = CM[i][j]=='0' ? ((CM[i+1][j]=='0'? 1:0)+(CM[i][j+1]=='0'? 1:0)+(CM[i][j-1]=='0'? 1:0)) : 0)==1)
        {
            dead_ends.emplace_back(i, j);
        }
        if(ACM[i][j]!=0)
        {
            open_cells++;
        }
    }
    
    for(int i=lines-1, j=1; j<columns-1; j++)
    {
        if((ACM[i][j] = CM[i][j]=='0' ? ((CM[i-1][j]=='0'? 1:0)+(CM[i][j+1]=='0'? 1:0)+(CM[i][j-1]=='0'? 1:0)) : 0)==1)
        {
            dead_ends.emplace_back(i, j);
        }
        if(ACM[i][j]!=0)
        {
            open_cells++;
        }
    }
    for(int i=1, j=0; i<lines-1; i++)
    {
        if((ACM[i][j] = CM[i][j]=='0' ? ((CM[i+1][j]=='0'? 1:0)+(CM[i-1][j]=='0'? 1:0)+(CM[i][j+1]=='0'? 1:0)) : 0)==1)
        {
            dead_ends.emplace_back(i, j);
        }
        if(ACM[i][j]!=0)
        {
            open_cells++;
        }
    }
    for(int i=1, j=columns-1; i<lines-1; i++)
    {
        if((ACM[i][j] = CM[i][j]=='0' ? ((CM[i+1][j]=='0'? 1:0)+(CM[i-1][j]=='0'? 1:0)+(CM[i][j-1]=='0'? 1:0)) : 0)==1)
        {
            dead_ends.emplace_back(i, j);
        }
        if(ACM[i][j]!=0)
        {
            open_cells++;
        }
    }

    //inner
    for(int i=1; i<lines-1; i++)
    {
        for(int j=1; j<columns-1; j++)
        {
            if((ACM[i][j] = CM[i][j]=='0' ? ((CM[i+1][j]=='0'? 1:0)+(CM[i-1][j]=='0'? 1:0)+(CM[i][j+1]=='0'? 1:0)+(CM[i][j-1]=='0'? 1:0)) : 0)==1)
            {
                dead_ends.emplace_back(i, j);
            }
            if(ACM[i][j]!=0)
            {
                open_cells++;
            }
        }
    }


    //dead-ends
    #pragma omp parallel for num_threads(max_threads)
    for(size_t i=0; i<dead_ends.size(); i++)
    {
        dead_end_thread(dead_ends[i]);
    }
    if(scide==2)
    {
        dece=true;
        sc = dead_end_sc(coords(start.line, start.column, lines, columns, start.line, start.column));
    }

    //ACM reduction: sets max visits to a cell based on its diagonals
    #pragma omp parallel for num_threads(max_threads)
    for(int i=1; i<lines-1; i++)
    {
        for(int j=1; j<columns-1; j++)
        {
            if(ACM[i][j]==4)
            {
                int t=0;
                if(CM[i+1][j+1]=='1')
                {
                    t++;
                }
                if(CM[i+1][j-1]=='1')
                {
                    t++;
                }
                if(CM[i-1][j+1]=='1')
                {
                    t++;
                }
                if(CM[i-1][j-1]=='1')
                {
                    t++;
                }
                if(t==0)
                {
                    ACM[i][j]=2;
                }
            }
        }
    }


    //add other procedures here

    //pathFinder
    thread p(progress);
    p.detach(); //progress saver

    walker(coords(start.line, start.column, lines, columns, start.line, start.column, sc, de_preppend));
    isrunning=false;
    pcv.notify_all();

    cout << "solution: ";
    for(size_t k=0; k<solution.size(); k++)
    {
        cout << solution[k];
    }
    cout << endl << "reg: " << endl;
    for(int k=0; k<lines; k++)
    {
        for(int l=0; l<columns; l++)
        {
            cout << s_reg[k][l];
        }
        cout << endl;
    }
    return 0;
}

void dead_end_thread(coords ic)
{
    int i=ic.line, j=ic.column;
    vector<int> path;
    while(true)
    {
        nav_mtxs[i][j].lock();
        if(i==start.line && j==start.column && ACM[i][j]<3)
        {
            scide=2;
        }
        else if(i==start.line && j==start.column && scide<=2)
        {
            scide+=1;
        }
        if(ACM[i][j]<3 && !(i==start.line && j==start.column))
        {
            int ti=i, tj=j;//values of i,j before updating
            CM[i][j]='1';
            if(dead_end_vector_position[i][j]!=0) //if some dead_end reaches this cell, append its path to this
            {
                path.insert(path.begin(), dead_ends_paths[dead_end_vector_position[i][j]-1].begin(), dead_ends_paths[dead_end_vector_position[i][j]-1].end());
            }
            if(i-1>=0 && CM[i-1][j]=='0')//up:0
            {
                i--;
                path.push_back(0);
                path.insert(path.begin(), 1);
                ACM[ti][tj]=0;
            }
            else if(i+1<lines && CM[i+1][j]=='0')//down:1
            {
                i++;
                path.push_back(1);
                path.insert(path.begin(), 0);
                ACM[ti][tj]=1;
            }
            else if(j-1>=0 && CM[i][j-1]=='0')//left:2
            {
                j--;
                path.push_back(2);
                path.insert(path.begin(), 3);
                ACM[ti][tj]=2;
            }
            else//right:3
            {
                j++;
                path.push_back(3);
                path.insert(path.begin(), 2);
                ACM[ti][tj]=3;
            }
            nav_mtxs[ti][tj].unlock();
        }
        else
        {
            ACM[i][j]--;
            if(dead_end_vector_position[i][j]!=0)
            {
                path.insert(path.begin(), dead_ends_paths[dead_end_vector_position[i][j]-1].begin(), dead_ends_paths[dead_end_vector_position[i][j]-1].end());
            }
            m_dep.lock();
            dead_ends_paths.push_back(path);
            dead_end_vector_position[i][j]=dead_ends_paths.size();
            m_dep.unlock();
            nav_mtxs[i][j].unlock();
            return;
        }
    }
}

int dead_end_sc(coords ic)
{
    int i=ic.line, j=ic.column;
    vector<int> spath, path; //spath: path with return for cell count calculations, path: path without return to be appended in the solution path
    while(true)
    {
        if(ACM[i][j]<3)
        {
            int ti=i, tj=j;//values of i,j before updating
            CM[i][j]='1';
            if(dead_end_vector_position[i][j]!=0) //if some dead_end reaches this cell, append its path to this
            {
                spath.insert(spath.begin(), dead_ends_paths[dead_end_vector_position[i][j]-1].begin(), dead_ends_paths[dead_end_vector_position[i][j]-1].end());
                path.insert(path.begin(), dead_ends_paths[dead_end_vector_position[i][j]-1].begin(), dead_ends_paths[dead_end_vector_position[i][j]-1].end());
            }
            if(i-1>=0 && CM[i-1][j]=='0')//up:0
            {
                i--;
                de_preppend.push_back(0);
                spath.push_back(0);
                spath.insert(spath.begin(), 1);
                path.insert(path.begin(), 1);
                ACM[ti][tj]=0;
            }
            else if(i+1<lines && CM[i+1][j]=='0')//down:1
            {
                i++;
                de_preppend.push_back(1);
                spath.push_back(1);
                spath.insert(spath.begin(), 0);
                path.insert(path.begin(), 0);
                ACM[ti][tj]=1;
            }
            else if(j-1>=0 && CM[i][j-1]=='0')//left:2
            {
                j--;
                de_preppend.push_back(2);
                spath.push_back(2);
                spath.insert(spath.begin(), 3);
                path.insert(path.begin(), 3);
                ACM[ti][tj]=2;
            }
            else//right:3
            {
                j++;
                de_preppend.push_back(3);
                spath.push_back(3);
                spath.insert(spath.begin(), 2);
                path.insert(path.begin(), 2);
                ACM[ti][tj]=3;
            }
        }
        else
        {
            ACM[i][j]--;
            start.line=i;
            start.column=j;
            de_append=path;
            return spath.size()/2;
        }
    }
}

void walker(coords ic)
{
    thread threads[4];
    int boolc;
    bool bools[4];
    
    if((solution.size()!=0 && ic.path.size()>=(size_t)solution.size()) || solution.size()==(size_t)open_cells) //BnB upper bound
    {
        return;
    }

    ic.reg[ic.line][ic.column]++;
    if(ic.direction!=-1)
    {
        ic.path.push_back(ic.direction);
    }

    if(dead_end_vector_position[ic.line][ic.column]!=0 && ic.reg[ic.line][ic.column]==1 && ic.path.size()!=de_preppend.size()) //append dead_ends
    {
        ic.path.insert(ic.path.end(), dead_ends_paths[dead_end_vector_position[ic.line][ic.column]-1].begin(), dead_ends_paths[dead_end_vector_position[ic.line][ic.column]-1].end());
        ic.cells+=dead_ends_paths[dead_end_vector_position[ic.line][ic.column]-1].size()/2;
    }

    if(ic.cells==open_cells && ic.line==start.line && ic.column==start.column) //caminho completo
    {
        m_dep.lock();
        if(solution.size()!=0 && ic.path.size()>=(size_t)solution.size()) //false positive
        {
            m_dep.unlock();
            return;
        }
        ic.path.insert(ic.path.end(), de_append.begin(), de_append.end());
        solution=ic.path;
        s_reg=ic.reg;
        m_dep.unlock();
        return;
    }

    bools[0]=ic.line-1>=0 && CM[ic.line-1][ic.column]=='0' && ic.reg[ic.line-1][ic.column]==0;
    bools[1]=ic.line+1<lines && CM[ic.line+1][ic.column]=='0' && ic.reg[ic.line+1][ic.column]==0;
    bools[2]=ic.column-1>=0 && CM[ic.line][ic.column-1]=='0' && ic.reg[ic.line][ic.column-1]==0;
    bools[3]=ic.column+1<columns && CM[ic.line][ic.column+1]=='0' && ic.reg[ic.line][ic.column+1]==0;
    boolc=bools[0]+bools[1]+bools[2]+bools[3]; //gives preference to never visited cells

    //calcular adjacencia
    int adj=(ic.line-1>=0 && ic.reg[ic.line-1][ic.column]==0)+(ic.line+1<lines && ic.reg[ic.line+1][ic.column]==0)+(ic.column-1>=0 && ic.reg[ic.line][ic.column-1]==0)+(ic.column+1>columns && ic.reg[ic.line][ic.column+1]==0);
    if(ic.reg[ic.line][ic.column]==1 && adj==0) //if it is the first time a cell is being visited, it has the right to return in the opposite direction
    {
        ic.direction=-1;
    }

    if(ic.line-1>=0 && CM[ic.line-1][ic.column]=='0' && ic.reg[ic.line-1][ic.column]<ACM[ic.line-1][ic.column] && ic.direction!=1 && (boolc==0||bools[0]))//up:0
    {
        if(runing_threads<max_threads)
        {
            try
            {
                runing_threads++;
                threads[0] = thread(thread(walker, coords(ic.line-1, ic.column, 0, ic.cells+(ic.reg[ic.line-1][ic.column]==0 ?1:0), ic.path, ic.reg)));
            }
            catch(exception& e)
            {
                runing_threads--;
                walker(coords(ic.line-1, ic.column, 0, ic.cells+(ic.reg[ic.line-1][ic.column]==0 ?1:0), ic.path, ic.reg));
            }
        }
        else
        {
            walker(coords(ic.line-1, ic.column, 0, ic.cells+(ic.reg[ic.line-1][ic.column]==0 ?1:0), ic.path, ic.reg));
        }
    }
    if(ic.line+1<lines && CM[ic.line+1][ic.column]=='0' && ic.reg[ic.line+1][ic.column]<ACM[ic.line+1][ic.column] && ic.direction!=0 && (boolc==0||bools[1]))//down:1
    {
        if(runing_threads<max_threads)
        {
            try
            {
                runing_threads++;
                threads[1] = thread(thread(walker, coords(ic.line+1, ic.column, 1, ic.cells+(ic.reg[ic.line+1][ic.column]==0 ?1:0), ic.path, ic.reg)));
            }
            catch(exception& e)
            {
                runing_threads--;
                walker(coords(ic.line+1, ic.column, 1, ic.cells+(ic.reg[ic.line+1][ic.column]==0 ?1:0), ic.path, ic.reg));
            }
        }
        else
        {
            walker(coords(ic.line+1, ic.column, 1, ic.cells+(ic.reg[ic.line+1][ic.column]==0 ?1:0), ic.path, ic.reg));
        }        
    }
    if(ic.column-1>=0 && CM[ic.line][ic.column-1]=='0' && ic.reg[ic.line][ic.column-1]<ACM[ic.line][ic.column-1] && ic.direction!=3 && (boolc==0||bools[2]))//left:2
    {
        if(runing_threads<max_threads)
        {
            try
            {
                runing_threads++;
                threads[2] = thread(thread(walker, coords(ic.line, ic.column-1, 2, ic.cells+(ic.reg[ic.line][ic.column-1]==0 ?1:0), ic.path, ic.reg)));
            }
            catch(exception& e)
            {
                runing_threads--;
                walker(coords(ic.line, ic.column-1, 2, ic.cells+(ic.reg[ic.line][ic.column-1]==0 ?1:0), ic.path, ic.reg));
            }
        }
        else
        {
            walker(coords(ic.line, ic.column-1, 2, ic.cells+(ic.reg[ic.line][ic.column-1]==0 ?1:0), ic.path, ic.reg));
        }        
    }
    if(ic.column+1<columns && CM[ic.line][ic.column+1]=='0' && ic.reg[ic.line][ic.column+1]<ACM[ic.line][ic.column+1] && ic.direction!=2 && (boolc==0||bools[3]))//right:3
    {
        if(runing_threads<max_threads)
        {
            try
            {
                runing_threads++;
                threads[3] = thread(thread(walker, coords(ic.line, ic.column+1, 3, ic.cells+(ic.reg[ic.line][ic.column+1]==0 ?1:0), ic.path, ic.reg)));
            }
            catch(exception& e)
            {
                runing_threads--;
                walker(coords(ic.line, ic.column+1, 3, ic.cells+(ic.reg[ic.line][ic.column+1]==0 ?1:0), ic.path, ic.reg));
            }
        }
        else
        {
            walker(coords(ic.line, ic.column+1, 3, ic.cells+(ic.reg[ic.line][ic.column+1]==0 ?1:0), ic.path, ic.reg));
        }        
    }

    for(size_t i=0; i<4; i++)//join threads
    {
        if(threads[i].joinable())
        {
            threads[i].join();
            runing_threads--;
        }
    }
    return;
}

void progress()
{
    while(isrunning)
    {
        m_dep.lock();
        pcv.wait_for(m_dep, chrono::seconds(60));
        if(!isrunning)
        {
            m_dep.unlock();
            return;
        }
        if(solution.size()==0)
        {
            m_dep.unlock();
            continue;
        }
        cout << "path: ";
        for(size_t k=0; k<solution.size(); k++)
        {
            cout << solution[k];
        }
        cout << endl << "reg: " << endl;
        for(int k=0; k<lines; k++)
        {
            for(int l=0; l<columns; l++)
            {
                cout << s_reg[k][l];
            }
            cout << endl;
        }
        cout << endl;
        m_dep.unlock();
    }
}
