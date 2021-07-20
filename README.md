# molecad

API для извлечения информации из баз данных Pubchem и дальнейшего взаимодействия с этими данными.

## Содержание

* [Документация на Confluence](#документация-на-confluence)
* [Настройка виртуального окружения](#настройка-виртуального-окружения)
* [Настройка переменных окружения](#настройка-переменных-окружения)
* [Использование консольной утилиты для скачивания свойств молекул из Pubchem](#использование-консольной-утилиты-для-скачивания-свойств-молекул-из-Pubchem)
* HOWTO
  * [Как развернуть MongoDB локально на macOS](#как-развернуть-mongodb-локально-на-macos)
    
* [Использование mongo-rdkit для подструктурного поиска по базе](#использование-mongo-rdkit-для-подструктурного-поиска-по-базе)
* [Схема базы данных](#схема-базы-данных)
* [Документация ручек API](#документация-ручек-api)
* [Структура проекта](#структура-проекта)


## Документация на Confluence

* [Как развернуть MongoDB локально на macOS](https://confluence.biocad.ru/x/3BhQCw)
* [Что можно скачать из баз данных PubChem](https://confluence.biocad.ru/x/DwpQCw)
* [Как использовать molecad для скачивания данных из Pubchem](https://confluence.biocad.ru/x/_ixQCw)
* [Описание ручек api](https://confluence.biocad.ru/x/jCxQCw)


## Настройка виртуального окружения

> Перед тем как приступить к настройке виртуального окружения, убедитесь, что на вашем компьютере 
> установлены [pyenv](https://github.com/pyenv/pyenv-installer) и [poetry](https://python-poetry.org/docs/), 
> причем их пути должны быть добавлены в `PATH`.

Склонируйте [репозиторий](https://github.com/vvrubel/molecad) и запустите shell из директории, 
содержащей файл `pyproject.toml`.

В проекте используется Python версии 3.7.10, который вы можете установить, используя `pyenv`:

```shell
$ pyenv install 3.7.10
$ pyenv local 3.7.10
```

Создать виртуальное окружение и установить все зависимости проекта можно с помощью `poetry`:

```shell
$ poetry install
```

Вуаля! Все зависимости проекта установлены и никакой конды, вам не потребовалось.  
Подробнее о реализации пути через `conda` можно почитать [тут](https://confluence.biocad.ru/x/vi1QCw).


## Настройка переменных окружения

После настройки виртуальной среды, необходимо задать переменные окружения – это параметры для подключения к mongodb и запуска скрипта, выполняющего запросы к базам данных Pubchem.  
Удобнее всего установить переменные один раз, создав файл `.env` – для тестовой базы данных и файл `prod.env`– для основной базы из шаблона `.env.sample`. 

В корневой директории найдите файл `.env.sample`, скопируйте его, заполните недостающие поля и переименуйте файл в `.env`; то же самое повторите для `prod.env`.

Переменные окружения могут быть двух типов: имеющие значение по умолчанию и не 
имеющие его. Если переменная имеет значение по умолчанию, то оно указано в таблице ниже. 
Переменные не имеющие значений по умолчанию обязательны к первоначальной настройке.

Переменная|Описание
---|---
MONGO_HOST|Хост базы данных MongoDB; по умолчанию -> `"127.0.0.1"`  
MONGO_PORT|Порт базы данных MongoDB; по умолчанию -> `27017`  
MONGO_USER|Имя пользователя базы данных MongoDB  
MONGO_PASSWORD|Пароль пользователя базы данных MongoDB  
MONGO_AUTH_SOURCE|База данных для аутентификации; по умолчанию -> `"admin"`  
MONGO_DB_NAME|Название рабочей базы данных; по умолчанию -> `"demo_db"`  
MONGO_PUBCHEM_COLLECTION|Название коллекции, в которой хранятся данные о свойствах молекул из базы данных `Compound`, скачанные с серверов Pubchem; по умолчанию -> `"pubchem"` 
MONGO_MOLECULES_COLLECTION|Название коллекции MongoDB, содержащая документы, являющиеся репрезентацией молекул – используется для подструктурного поиска; по умолчанию -> `"molecules"`  
MONGO_MFP_COUNTS_COLLECTION|Название коллекции MongoDB, содержащая документы, являющиеся репрезентацией молекул – используется для подструктурного поиска; по умолчанию -> `"mfp_counts"`  
MONGO_PERMUTATIONS_COLLECTION|Название коллекции MongoDB, содержащая документы, являющиеся репрезентацией молекул – используется для подструктурного поиска; по умолчанию -> `"permutations"`  
PROJ_DIR|Путь до корневой директории проекта, если значение не определено, то будет использоваться текущая директория; по умолчанию -> `.`
FETCH_DIR|Шаблон именования пути до директорий для сохранения файлов, полученных от серверов Pubchem с помощью команды `fetch`; по умолчанию -> `./data/fetch`
SPLIT_DIR|Шаблон именования пути до директорий для сохранения файлов, полученных в результате выполнения команды `split` над json-файлами; по умолчанию -> `./data/split`

## Использование консольной утилиты для скачивания свойств молекул из Pubchem

Данный проект предоставляет доступ к консольной утилите для извлечения информации о свойствах молекул с серверов Pubchem, сохранения полученных данных в файлы формата JSON, а также загрузки содержимого в базу данных MongoDB. 

Справка по работе утилиты:

```shell
$ poetry run python -m molecad.cli --help  

Usage: python -m molecad.cli [OPTIONS] COMMAND [ARGS]...

  Консольная утилита для извлечения информации о свойствах молекул с серверов
  Pubchem, сохранения полученной информации в файлы формата JSON. Также с
  помощью утилиты можно загружать полученные данные в базу MongoDB и
  подготавливать их для поиска молекул.Для выполнения команд, которые
  устанавливают соединение с базой данных, необходимо указать опцию
  --prod/--dev для выбора рабочей базы.

Options:
  --setup [PROD|DEV]  Опция, позволяющая выбирать конфигурационный файл,
                      содержащий переменные окружения, который также
                      определяет настройки базы данных.
  --help              Show this message and exit.


Commands:
  fetch     Извлекает и сохраняет информацию из базы данных Pubchem –...
  populate  Из указанной директории загружает файлы в коллекцию MongoDB,...
  split     При использовании команды `db.collection.insert_many({..})`...

```

При скачивании данных с серверов Pubchem с помощью команды `fetch` консольной утилиты, 
программа не ждет пока выполнятся все запросы для того, чтобы сохранить все данные в файл, 
а пишет результаты по мере выполнения.  

В ходе работы команды `fetch` создается поддиректория, называемая в соответствие с 
параметрами `start` и `stop`, в директории `./data/fetch/`; в созданную поддиректорию 
сохраняются файлы, имеющие названия, составленные из крайних значений диапазонов идентификаторов 
сохраняемого файла. Размер диапазона задается при вызове команды `fetch` с опцией `--f-size`. 

Так как MongoDB имеет ограничение на количество создаваемых документов при разовой загрузке в базу 
из файла, то дальнейший план будет определяться длинной списка, который находится в этом файле. 
Если длина списка (число запрошенных идентификаторов) < 100000, то этот файл можно сразу загрузить 
в локальную базу данных MongoDB, используя команду `populate`. Иначе перед загрузкой потребуется 
его разделение на chunked-файлы меньшего размера с помощью команды `split`, например по 1000 
идентификаторов в файле.

Запуск команд утилиты рекомендуется производить из корневой директории проекта. Команды, которые 
непосредственно взаимодействуют с базой данных, требуют наличие запущенной службы `mongod`. Если 
вы уже развернули локальную базу данных, то запустить сервер можно с помощью команды:
```shell
mongod [--dbpath <path>] --auth
```
где <path> – путь до локальной базы MongoDB.


Справка по команде `fetch`:  
```shell
$ poetry run python -m molecad.cli --setup DEV fetch --help
  
Выбрана база demo_db
Usage: python -m molecad.cli fetch [OPTIONS]

  Извлекает и сохраняет информацию из базы данных Pubchem – 'Compound'. Для
  совершения запроса к серверу необходимо уточнить диапазон идентификаторов
  интересуемых соединений - 'CID'. Данные извлекаются путем отправления
  запроса на сгенерированный URL. Из-за ограничения на количество символов в
  строке URL, отправка запроса, включающего все желаемые идентификаторы, может
  привести к ошибке, поэтому рекомендуется использовать опцию `--size`,
  которая по умолчанию равна `100`. Полученные данные сохраняются в файл,
  причем имеется возможность сохранять файл порциями до окончания работы
  программы, указав в качестве опции `--f-size` желаемое количество
  идентификаторов в сохраняемом файле, которое по умолчанию равна `1000`.
  Пример запуска команды:
  $ poetry run python -m molecad.cli --setup DEV fetch 
                                     --start 1
                                     --stop 1000000 
                                     --size 100 
                                     --f-size 1000  

Options:
  --start TEXT      Первое значение из запрашиваемых CID.  [required]
  --stop TEXT       Последнее значение из запрашиваемых CID.  [required]
  --size INTEGER    Максимальное число идентификаторов в одном запросе, по
                    умолчанию равно 100.  [default: 100]
  --f-size INTEGER  Максимальное число идентификаторов в сохраняемом файле, по
                    умолчанию равно 1000.  [default: 1000]
  --help            Show this message and exit.
```

Справка по команде `split`:

```shell
$ poetry run python -m molecad.data.console split --help

Usage: python -m molecad.data.console.cli split [OPTIONS]

  При использовании команды `db.collection.insert_many({..})` имеется
  ограничение на максимально допустимое количество добавляемых документов за
  один раз равное 100000. Данная функция служит для того, чтобы разрезать
  JSON-файлы, превышающие указанное выше ограничение, на файлы меньшего
  размера.
  Пример запуска команды:
  $ poetry run python -m molecad.cli --setup DEV split --file path/to/file.json \
                                                       --f-size 1000


Options:
  --file PATH       Файл не подходящий под критерии загрузки файлов в MongoDB.
                    [required]
  --f-size INTEGER  Максимальное число идентификаторов в сохраняемом файле, по
                    умолчанию равно 1000.  [default: 1000]
  --help            Show this message and exit.
```

Справка по команде `populate`, обратите внимание, что необходимо указать опцию `--setup` для выбора настроек `PROD` или `DEV` окружения:

```shell
$ poetry run python -m molecad.cli --setup DEV populate --help
           
Выбрана база demo_db
Usage: python -m molecad.cli populate [OPTIONS]

  Из указанной директории загружает файлы в коллекцию MongoDB, при условии что
  файл содержит не более 100000 документов. Перед загрузкой документов на
  коллекции будут созданы уникальные индекс 'CID'.
  Пример запуска команды:
  poetry run python -m molecad.cli 
             --setup DEV populate 
             --f-dir ./data/fetch
             --drop

Options:
  --f-dir PATH  Путь до директории, содержащей chunked-файлы, каждый из
                которых представляет собой список длиной до 100000 элементов.
                [required]
  --drop        Если опция '--drop' указана, то коллекция будет очищена.
  --help        Show this message and exit.
```

## Использование mongo-rdkit для подструктурного поиска по базе

Исторически сложилось, что `mongo-rdkit` и просто `rdkit` используют виртуальное окружение 
`conda` для разрешения своих зависимостей.  Но в данном проекте было решено реализовать 
управление зависимостями этих библиотек с использованием `poetry`.  
Если вы хотите использовать окружение конды и вам хочется узнать, как установить `conda` рядом с 
`pyenv` наиболее безболезненно – читайте [тут](https://confluence.biocad.ru/x/vi1QCw).

После того как мы настроили окружение, его переменные, развернули MongoDB и загрузили все 
скачанные данные в локальную базу – можно приступать к рассмотрению задачи о подструктурном 
поиске. Перед этим нужно отметить несколько ключевых моментов, которые были учтены при 
выполнении команды `populate`:

* По причине того что поиск производится по smiles молекул, каждый документ должен 
  содержать поле `CanonicalSMILES`. В противном случае такие документы будут удалены из коллекции.
* При загрузке в локальную базу данных, извлеченных из Pubchem, первичный индекс `_id` документов 
  генерируется самостоятельно. Это может повлечь за собой то, что в случае повторной загрузки 
  документов база сгенерирует для них новый уникальный `_id`, не сообщив об этом пользователю. 
  Во избежание подобных ситуаций – дублирования документов – необходимо перед загрузкой документов 
  создать на коллекции новый уникальный индекс, который привязан к уникальному полю документа. 
  В нашем случае им является поле `CID`.
* `rdkit` или `mongo-rdkit` не генерирует схемы для металлогранических соединений, поэтому мы 
  тоже удаляем такие документы.



Модуль `Database.registration` создает документальное представление молекул в соответствии с 
генерируемой схемой и обрабатывает настройки их записи в базу данных. 
Этот модуль состоит из двух частей. 
Во-первых, он определяет глобальную переменную `HASH_FUNCTIONS` как словарь, который отображает 
имена хэш-функций на методы. Он также определяет глобальные переменные `DEFAULT_SCHEME_NAME`, 
`DEFAULT_AUTHOR`, `DEFAULT_PREPROCESS` и `DEFAULT_INDEX`, которые используются при создании схемы и, 
таким образом, определены для упрощения настройки. 
Во-вторых, файл определяет объект `MolDocScheme`, который содержит информацию о схеме переменных 
экземпляра и передается в методы `.write` для создания документа молекулы. 
По умолчанию `MolDocScheme` включает имя схемы, автора, независимо от того, была ли молекула 
предварительно обработана, параметр индекса, два хэша, отпечатки пальцев и поля значений. 
Вся информация, содержащаяся в объекте `MolDocScheme`, может использоваться непосредственно для 
создания документов для молекул:
```python
rdmol = Chem.MolFromSmiles('Cc1ccccc1')
scheme = registration.MolDocScheme()
scheme.generate_mol_doc(rdmol)
```
Класс `MolDocScheme` также определяет ряд методов экземпляра, таких как `MolDocScheme.set_index` и 
`MolDocScheme.remove_field`, которые можно использовать для изменения схем документов:
```python
scheme.remove_field('CanonicalSmiles')
scheme.add_hash_field('MolFormula')
scheme.set_index('MolFormula')
scheme.generate_mol_doc(rdmol)
```
Поскольку объекты `MolDocScheme` не содержат функций, а только ссылки на них – их можно обработать. 
Фактически, методы записи могут сохранять `MolDocScheme`s, таким образом схемы, созданные 
пользователем ранее, могут быть в последующем также извлечены. 
В случае, если пользователи не работают с SDF, `.write` также предоставляет `WriteFromMolList`, 
который будет принимать список Python объектов `rdmol` вместо аргумента `SDF` в `WriteFromSDF`.


## HOWTO

### Как развернуть MongoDB локально на macOS

Установите на свой компьютер сервер MongoDB: Community Server можно поставить с помощью 
менеджера пакетов для macOS - `brew`, а Enterprise Server - из архива `.targz`. Добавьте путь до 
исполняемого файла в `PATH`.  

Ниже описаны основные шаги разворачивания локальной базы на примере Enterprise Server. Детальное 
описание процесса установки и инструкции для Community Server можно найти по [ссылке](https://confluence.biocad.ru/x/3BhQCw).

Переходим в каталог с только что установленным сервером MongoDB и создаем в нем подкаталог для 
локальной базы, после чего запустим сервер.

```shell
$ mkdir -p ./data/db
$ mongod --dbpath ./data/db
```

После того как сервер запущен, подключаемся к локальной базе данных. При первом запуске базы необходимо создать пользователей.  

Пример кода для создания суперпользователя: 

```shell
$ mongo
MongoDB Enterprise> use admin
MongoDB Enterprise> db.createUser({ user: "superuser", pwd: "<pwd>", roles: [ "root" ] })
```

После того как пользователь будет создан, нужно выключить сервер и выйти в исходный shell.

```shell
MongoDB Enterprise> db.shutdownServer()
MongoDB Enterprise> exit
```

Теперь перезапускаем сервер mongod, но с параметром `--auth` для подключения по логину-паролю.
Готово! Теперь вы можете использовать свою базу.

```shell
$ mongod --dbpath ./data/db --auth
$ mongo mongodb://superuser:<password>@127.0.0.1:27017/admin
```
#### Создание конфигурационного файла

```shell
mongod --config /etc/mongod.conf
mongod -f /etc/mongod.conf
```

#### Как проверить запущен ли сервер mongod

```shell
ps aux | grep -v grep | grep mongod
```

## Схема базы данных
```
    CID: int
    MolecularFormula: str
    MolecularWeight: Union[float, str]  # почему-то в response приходит строка
    CanonicalSMILES: str
    InChI: str
    IUPACName: str
    XLogP: float
    HBondDonorCount: int
    HBondAcceptorCount: int
    RotatableBondCount: int
    AtomStereoCount: int
    BondStereoCount: int
    Volume3D: Union[float, int]
```
## Документация ручек API