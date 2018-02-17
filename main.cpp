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

#define DEFAULT_MAX_THREADS 2

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

    coords(int a=0, int b=0, int l=0, int c=0, int sl=0, int sc=0)
    {
        line=a;
        column=b;
        direction=-1;
        cells=0;
        reg=vector<vector<int>>(l, vector<int>(c));
        if(l!=0 && c!=0)
        {
            reg[sl][sc]=-1;
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

class scheduler
{
    private:
    queue<void (*)(coords)> functions;
    queue<coords> parameters;
    int runing_threads, scheduled_threads, dispatched_threads;
    mutex m;

    template <class T> void clear( queue<T> &q )
    {
    std::queue<T> empty;
    std::swap( q, empty );
    }

    public:
    int max_threads;
    condition_variable_any nrt; 
    mutex m_mtx;

    scheduler(int mt)
    {
        max_threads=mt;
        runing_threads=0;
        scheduled_threads=0;
        dispatched_threads=0;
    }

    ~scheduler()
    {
        clear<void (*)(coords)>(functions);
        clear<coords>(parameters);
    }

    void set_scheduled_threads(int a)
    {
        scheduled_threads=a;
        dispatched_threads=0;
    }

    void overwrite_scheduled_threads(int a)
    {
        m.lock();
        scheduled_threads+=a;
        m.unlock();
    }

    void schedule(void (*f)(coords), coords c)
    {
        m.lock();
        if(runing_threads<max_threads)
        {
            runing_threads++;
            thread t(f, c);
            t.detach();
        }
        else
        {
            functions.push(f);
            parameters.push(c);            
        }
        m.unlock();
    }

    void dispatch()
    {
        m.lock();
        runing_threads--;
        dispatched_threads++;
        if(parameters.size()>0)
        {
            runing_threads++;
            thread t(functions.front(), parameters.front());
            t.detach();
            functions.pop();
            parameters.pop();            
        }
        if(threads_completed()) nrt.notify_one();
        m.unlock();
    }

    bool threads_completed()
    {
        return dispatched_threads==scheduled_threads ? true : false;
    }

    bool free_threads()
    {
        return (runing_threads<max_threads);
    }

    void print_thread_info()
    {
        cout << "st: " << scheduled_threads << " rt: " << runing_threads << " dp: " << dispatched_threads << " queued_t: " << scheduled_threads-dispatched_threads << endl;
    }
};

void dead_end_thread(coords);
void walker(coords);

scheduler* S;

vector<string> CM;
vector<vector<int>> ACM;
vector<vector<int>> BnB;
vector<vector<int>> dead_end_vector_position; //fixed(can be improved)
vector<vector<int>> dead_ends_paths;//fixed(can be improved)
coords start;
vector<int> solution;//decreasing
int lines=0, columns=0, open_cells=0, max_threads, runing_threads=0;

mutex** nav_mtxs;//fixed
mutex m_dep;

int main(int argc, char* argv[])
{
    S=new scheduler(argc>1 ? stoi(argv[1]) : DEFAULT_MAX_THREADS);
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


    //conversion
    ACM.resize(lines);
    BnB.resize(lines);
    dead_end_vector_position.resize(lines);
    nav_mtxs  = (mutex**)malloc(sizeof(mutex*)*lines);
    columns=CM[0].size();
    for(int i=0; i<lines; i++)
    {
        dead_end_vector_position[i].resize(columns);
        nav_mtxs [i] = (mutex*)malloc(sizeof(mutex)*columns);
        ACM[i].resize(columns);
        BnB[i].resize(columns);
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
    S->set_scheduled_threads(dead_ends.size());
    for(size_t i=0; i<dead_ends.size(); i++)
    {
        S->schedule(dead_end_thread, dead_ends[i]);
    }
    if(dead_ends.size()>0)
    {
        S->m_mtx.lock();//avoid undefined behavior
        S->nrt.wait(S->m_mtx);
        S->m_mtx.unlock();
    }

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
    //if ic is inside dead_end, follow ACM to find the exit TODO:v2

    //pathFinder
    walker(coords(start.line, start.column, lines, columns, start.line, start.column));
    //cout << "results have been saved to [arquivo]. \nPress any key to continue...";
    cout << "solution: ";
    for(size_t k=0; k<solution.size(); k++)
    {
        cout << solution[k];
    }
    cout << endl;
}

void dead_end_thread(coords ic)
{
    int i=ic.line, j=ic.column;
    vector<int> path;
    while(true)
    {
        nav_mtxs[i][j].lock();
        if(ACM[i][j]<3)
        {
            int ti=i, tj=j;//values of i,j before updating
            CM[i][j]='1';
            if(dead_end_vector_position[i][j]!=0)
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
            S->dispatch();
            return;
        }
    }
}

void walker(coords ic)
{
    thread threads[4];
    
    if((solution.size()!=0 && ic.path.size()>=(size_t)solution.size()) || solution.size()==(size_t)open_cells) //BnB upper bound
    {
        return;
    }

    ic.reg[ic.line][ic.column]++;
    if(ic.direction!=-1)
    {
        ic.path.push_back(ic.direction);
    }

    if(dead_end_vector_position[ic.line][ic.column]!=0 && ic.reg[ic.line][ic.column]<=1) //append dead_ends
    {
        ic.path.insert(ic.path.end(), dead_ends_paths[dead_end_vector_position[ic.line][ic.column]-1].begin(), dead_ends_paths[dead_end_vector_position[ic.line][ic.column]-1].end());
        ic.cells+=dead_ends_paths[dead_end_vector_position[ic.line][ic.column]-1].size()/2;
    }

    if(ic.cells==open_cells && ic.line==start.line && ic.column==start.column) //caminho completo
    {
        m_dep.lock();
        if(solution.size()!=0 && ic.path.size()>=(size_t)solution.size())
        {
            m_dep.unlock();
            return;
        }
        /*cout << "size: " << ic.path.size() << " path: ";
        for(size_t k=0; k<ic.path.size(); k++)
        {
            cout << ic.path[k];
        }
        cout << endl << "reg:\n";
        for(int k=0; k<lines; k++)
        {
            for(int l=0; l<columns; l++)
            {
                cout << ic.reg[k][l];
            }
            cout << endl;
        }
        cout << endl;*/
        solution=ic.path;
        m_dep.unlock();
        return;
    }

    if(ic.line-1>=0 && CM[ic.line-1][ic.column]=='0' && ic.reg[ic.line-1][ic.column]<ACM[ic.line-1][ic.column] && ic.direction!=1)//up:0
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
    if(ic.line+1<lines && CM[ic.line+1][ic.column]=='0' && ic.reg[ic.line+1][ic.column]<ACM[ic.line+1][ic.column] && ic.direction!=0)//down:1
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
    if(ic.column-1>=0 && CM[ic.line][ic.column-1]=='0' && ic.reg[ic.line][ic.column-1]<ACM[ic.line][ic.column-1] && ic.direction!=3)//left:2
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
    if(ic.column+1<columns && CM[ic.line][ic.column+1]=='0' && ic.reg[ic.line][ic.column+1]<ACM[ic.line][ic.column+1] && ic.direction!=2)//right:3
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
    for(size_t i=0; i<4; i++)
    {
        if(threads[i].joinable())
        {
            threads[i].join();
            runing_threads--;
        }
    }
    return;
}