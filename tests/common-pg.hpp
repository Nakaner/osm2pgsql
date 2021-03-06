#ifndef OSM2PGSQL_TEST_COMMON_PG_HPP
#define OSM2PGSQL_TEST_COMMON_PG_HPP

#include <cstdio>
#include <stdexcept>
#include <string>

#include <libpq-fe.h>

#include <boost/lexical_cast.hpp>

#include "format.hpp"
#include "options.hpp"
#include <catch.hpp>

#ifdef _MSC_VER
#include <process.h>
#include <windows.h>
#define getpid _getpid
#endif

/// Helper classes for postgres connections
namespace pg {

class result_t
{
public:
    result_t(PGresult *result) : m_result(result) {}

    ~result_t() noexcept { PQclear(m_result); }

    ExecStatusType status() const noexcept { return PQresultStatus(m_result); }

    int num_tuples() const noexcept { return PQntuples(m_result); }

    std::string get_value(int row, int col) const noexcept
    {
        return PQgetvalue(m_result, row, col);
    }

    bool is_null(int row, int col) const noexcept
    {
        return PQgetisnull(m_result, row, col) != 0;
    }

private:
    PGresult *m_result;
};

class conn_t
{
public:
    conn_t(char const *conninfo)
    {
        m_conn = PQconnectdb(conninfo);

        if (PQstatus(m_conn) != CONNECTION_OK) {
            fprintf(stderr, "Could not connect to database '%s' because: %s\n",
                    conninfo, PQerrorMessage(m_conn));
            PQfinish(m_conn);
            throw std::runtime_error("Database connection failed");
        }
    }

    ~conn_t() noexcept
    {
        if (m_conn) {
            PQfinish(m_conn);
        }
    }

    void exec(std::string const &cmd, ExecStatusType expect = PGRES_COMMAND_OK)
    {
        result_t res = query(cmd);
        if (res.status() != expect) {
            fprintf(stderr, "Query '%s' failed with: %s\n", cmd.c_str(),
                    PQerrorMessage(m_conn));
            throw std::runtime_error("DB exec failed.");
        }
    }

    result_t query(std::string const &cmd) const
    {
        return PQexec(m_conn, cmd.c_str());
    }

    template <typename T>
    T require_scalar(std::string const &cmd) const
    {
        result_t res = query(cmd);
        REQUIRE(res.status() == PGRES_TUPLES_OK);
        REQUIRE(res.num_tuples() == 1);

        std::string str = res.get_value(0, 0);
        return boost::lexical_cast<T>(str);
    }

    void assert_double(double expected, std::string const &cmd) const
    {
        REQUIRE(Approx(expected).epsilon(0.01) == require_scalar<double>(cmd));
    }

    void assert_null(std::string const &cmd) const
    {
        result_t res = query(cmd);
        REQUIRE(res.status() == PGRES_TUPLES_OK);
        REQUIRE(res.num_tuples() == 1);
        REQUIRE(res.is_null(0, 0));
    }

    result_t require_row(std::string const &cmd) const
    {
        result_t res = query(cmd);
        REQUIRE(res.status() == PGRES_TUPLES_OK);
        REQUIRE(res.num_tuples() == 1);

        return res;
    }

    unsigned long get_count(char const *table_name,
                            std::string const &where = "") const
    {
        auto const query = "SELECT count(*) FROM {} {} {}"_format(
            table_name, (where.empty() ? "" : "WHERE"), where);

        return require_scalar<unsigned long>(query);
    }

    void require_has_table(char const *table_name) const
    {
        auto const where = "oid = '{}'::regclass"_format(table_name);

        REQUIRE(get_count("pg_catalog.pg_class", where) == 1);
    }

private:
    PGconn *m_conn = nullptr;
};

class tempdb_t
{
public:
    tempdb_t()
    {
        try {
            conn_t conn("dbname=postgres");

            m_db_name = "osm2pgsql-test-{}-{}"_format(getpid(), time(nullptr));
            conn.exec("DROP DATABASE IF EXISTS \"{}\""_format(m_db_name));
            conn.exec("CREATE DATABASE \"{}\" WITH ENCODING 'UTF8'"_format(
                m_db_name));

            conn_t local = connect();
            local.exec("CREATE EXTENSION postgis");
            local.exec("CREATE EXTENSION hstore");
        } catch (std::runtime_error const &e) {
            fprintf(stderr, "Test database cannot be created: %s\n", e.what());
            fprintf(stderr, "Did you mean to run 'pg_virtualenv ctest'?\n");
            exit(1);
        }
    }

    ~tempdb_t() noexcept
    {
        if (!m_db_name.empty()) {
            try {
                conn_t conn("dbname=postgres");
                conn.query("DROP DATABASE IF EXISTS \"{}\""_format(m_db_name));
            } catch (...) {
                fprintf(stderr, "DROP DATABASE failed. Ignored.\n");
            }
        }
    }

    conn_t connect() const { return conn_t(conninfo().c_str()); }

    std::string conninfo() const { return "dbname=" + m_db_name; }

    database_options_t db_options() const
    {
        database_options_t opt;
        opt.db = m_db_name;

        return opt;
    }

private:
    std::string m_db_name;
};

} // namespace pg

#endif // OSM2PGSQL_TEST_COMMON_PG_HPP
